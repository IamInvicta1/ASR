#!/usr/bin/env python3
import os
import sys
import time
import csv
import subprocess
from typing import List, Tuple, Dict
import numpy as np


# -----------------------------
# CNF parsing / writing
# -----------------------------
def read_dimacs(path: str) -> Tuple[int, List[List[int]], List[str]]:
    """
    Returns (nvars, clauses, comments)
    clauses: list of list of ints (literals)
    """
    comments = []
    nvars = 0
    nclauses = 0
    clauses: List[List[int]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("c"):
                comments.append(line)
                continue
            if line.startswith("p"):
                parts = line.split()
                if len(parts) >= 4 and parts[1] == "cnf":
                    nvars = int(parts[2])
                    nclauses = int(parts[3])
                continue
            # clause line
            parts = line.split()
            lits = [int(x) for x in parts]
            if not lits:
                continue
            if lits[-1] == 0:
                lits = lits[:-1]
            # allow empty clause (edge case)
            clauses.append(lits)
    # nclauses in header might differ; we trust actual clauses list length
    if nvars <= 0:
        raise ValueError("No valid DIMACS header found (p cnf ...).")
    return nvars, clauses, comments


def write_dimacs(path: str, nvars: int, clauses: List[List[int]], comments: List[str] = None) -> None:
    if comments is None:
        comments = []
    with open(path, "w") as f:
        for c in comments:
            f.write(c + "\n")
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cl in clauses:
            if len(cl) == 0:
                f.write("0\n")
            else:
                f.write(" ".join(str(x) for x in cl) + " 0\n")


def add_unit_clauses(src_cnf: str, dst_cnf: str, new_units: List[int]) -> None:
    nvars, clauses, comments = read_dimacs(src_cnf)
    for lit in new_units:
        clauses.append([int(lit)])
    write_dimacs(dst_cnf, nvars, clauses, comments)


# -----------------------------
# Glucose runner + log parsing
# -----------------------------
def run_glucose(glucose_bin: str, cnf_path: str, log_path: str) -> None:
    # capture stdout+stderr to log file
    with open(log_path, "w") as f:
        subprocess.run([glucose_bin, cnf_path], stdout=f, stderr=subprocess.STDOUT, check=False)


def parse_glucose_log(log_path: str) -> Tuple[str, int, float]:
    """
    Returns (status, conflicts, cpu_time)
    status in {"SAT","UNSAT","UNKNOWN"}
    """
    status = "UNKNOWN"
    conflicts = None
    cpu_time = None

    with open(log_path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("s "):
                if "SATISFIABLE" in s and "UNSAT" not in s:
                    status = "SAT"
                elif "UNSATISFIABLE" in s:
                    status = "UNSAT"
            # conflicts line in Glucose typically: "c conflicts             : 206 ..."
            if s.lower().startswith("c conflicts"):
                # split on ":" then take first number
                parts = s.split(":")
                if len(parts) >= 2:
                    rhs = parts[1].strip().split()
                    if rhs:
                        try:
                            conflicts = int(rhs[0])
                        except:
                            pass
            if s.lower().startswith("c cpu time"):
                parts = s.split(":")
                if len(parts) >= 2:
                    rhs = parts[1].strip().split()
                    if rhs:
                        try:
                            cpu_time = float(rhs[0])
                        except:
                            pass

    if conflicts is None:
        conflicts = -1
    if cpu_time is None:
        cpu_time = -1.0
    return status, conflicts, cpu_time


# -----------------------------
# Rayleigh-proxy picker (v4)
# -----------------------------
def build_incidence(clauses: List[List[int]], nvars: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build incidence lists:
      clause_vars[j] = list of variable indices (0-based) in clause j
      var_clauses[i] = list of clause indices where var i appears (either sign)
    """
    clause_vars: List[List[int]] = []
    var_clauses: List[List[int]] = [[] for _ in range(nvars)]
    for j, cl in enumerate(clauses):
        vs = []
        for lit in cl:
            v = abs(lit) - 1
            if 0 <= v < nvars:
                vs.append(v)
                var_clauses[v].append(j)
        clause_vars.append(vs)
    return np.array(clause_vars, dtype=object), np.array(var_clauses, dtype=object)


def power_iter_variable_scores(clause_vars: np.ndarray, nvars: int, iters: int = 20) -> np.ndarray:
    """
    Compute a cheap global variable centrality score via power iteration on
    co-occurrence graph G = B^T B (implicit), where B is clause-variable incidence.

    We never build G explicitly:
      x_{t+1}[i] = sum_{clauses containing i} sum_{vars in that clause} x_t[var]
    """
    x = np.ones(nvars, dtype=np.float64)
    x /= np.linalg.norm(x) + 1e-12

    for _ in range(iters):
        y = np.zeros(nvars, dtype=np.float64)
        # for each clause, accumulate its mass to all vars in clause
        for vs in clause_vars:
            if len(vs) == 0:
                continue
            s = float(np.sum(x[vs]))
            y[vs] += s
        norm = np.linalg.norm(y) + 1e-12
        x = y / norm
    return x


def literal_rayleigh_proxy_scores(clauses: List[List[int]], nvars: int) -> Dict[int, float]:
    """
    Compute a Rayleigh-style proxy benefit per literal.
    We:
      1) compute global var scores x by power iteration (global coherence)
      2) compute clause weight w_j = sum_{v in clause j} x[v]
      3) define benefit(lit) ~ sum_{clauses satisfied by lit} w_j^2  - alpha * sum_{clauses where ~lit appears} w_j
         (removing satisfied clauses strongly reduces structure; falsifying just removes one edge)

    Returns dict: lit -> score (higher is better)
    """
    # parse cnf
    clause_vars, _ = build_incidence(clauses, nvars)
    x = power_iter_variable_scores(clause_vars, nvars, iters=25)

    # clause weight
    w = np.zeros(len(clauses), dtype=np.float64)
    for j, vs in enumerate(clause_vars):
        if len(vs) == 0:
            w[j] = 0.0
        else:
            w[j] = float(np.sum(x[vs]))

    alpha = 0.15  # small penalty factor

    score: Dict[int, float] = {}
    # initialize all literals
    for v in range(1, nvars + 1):
        score[v] = 0.0
        score[-v] = 0.0

    for j, cl in enumerate(clauses):
        wj = float(w[j])
        wj2 = wj * wj
        for lit in cl:
            # satisfied if we set lit True -> remove this clause
            score[lit] += wj2
            # falsifying opposite literal removes just one edge: small penalty
            score[-lit] -= alpha * wj

    return score


def pick_units_rayleigh(cnf_path: str, k: int) -> List[int]:
    nvars, clauses, _ = read_dimacs(cnf_path)
    sc = literal_rayleigh_proxy_scores(clauses, nvars)

    # Choose best literal per variable (avoid duplicates variable both signs)
    best_per_var: List[Tuple[float, int]] = []
    for v in range(1, nvars + 1):
        sp = sc.get(v, 0.0)
        sn = sc.get(-v, 0.0)
        if sp >= sn:
            best_per_var.append((sp, v))
        else:
            best_per_var.append((sn, -v))

    best_per_var.sort(reverse=True, key=lambda t: t[0])
    chosen = [lit for _, lit in best_per_var[:k]]
    return chosen


# -----------------------------
# Main iterative runner
# -----------------------------
def ensure_dirs():
    os.makedirs("shared/cnf", exist_ok=True)
    os.makedirs("shared/logs", exist_ok=True)
    os.makedirs("shared/results", exist_ok=True)


def main():
    if len(sys.argv) != 6:
        print("Usage:")
        print("  python3 versions/asr_v4_rayleigh/run_asr_v4_rayleigh.py INPUT_CNF TAG MAX_ITERS K_PER_ITER GLUCOSE_BIN")
        print("Example:")
        print("  python3 versions/asr_v4_rayleigh/run_asr_v4_rayleigh.py shared/cnf/n120_s2.cnf n120_s2_v4 30 1 ~/glucose/simp/glucose")
        sys.exit(1)

    inp = sys.argv[1]
    tag = sys.argv[2]
    max_iters = int(sys.argv[3])
    k_per_iter = int(sys.argv[4])
    glucose_bin = sys.argv[5]

    ensure_dirs()

    csv_path = f"shared/results/{tag}_asr_v4_rayleigh.csv"
    curr_cnf = inp

    # CSV header
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iter", "k_added_total", "new_units", "status", "conflicts", "cpu_time", "cnf_path", "log_path"])

    k_total = 0

    for it in range(max_iters + 1):
        cnf_out = f"shared/cnf/{tag}_iter{it}.cnf"
        log_out = f"shared/logs/{tag}_iter{it}.log"

        if it == 0:
            # just copy original cnf into iter0 for traceability
            nvars, clauses, comments = read_dimacs(curr_cnf)
            write_dimacs(cnf_out, nvars, clauses, comments)
        else:
            new_units = pick_units_rayleigh(curr_cnf, k_per_iter)
            k_total += len(new_units)
            add_unit_clauses(curr_cnf, cnf_out, new_units)

        # run solver
        run_glucose(glucose_bin, cnf_out, log_out)
        status, conflicts, cpu_time = parse_glucose_log(log_out)

        # record
        if it == 0:
            new_units_str = ""
        else:
            new_units_str = " ".join(str(x) for x in new_units)

        with open(csv_path, "a", newline="") as f:
            w = csv.writer(f)
            w.writerow([it, k_total if it > 0 else 0, new_units_str, status, conflicts, cpu_time, cnf_out, log_out])

        # terminal print
        if it == 0:
            print(f"[iter {it}] status={status} conflicts={conflicts} cpu={cpu_time}")
        else:
            print(f"[iter {it}] +{len(new_units)} units (total={k_total}): [{new_units_str}] | status={status} conflicts={conflicts} cpu={cpu_time}")

        # stop conditions
        if conflicts == 0 and status in ("SAT", "UNSAT"):
            print("Solved with 0 conflicts. Stopping.")
            break
        if status in ("SAT", "UNSAT") and conflicts >= 0 and conflicts < 5:
            # optional early exit when essentially collapsed
            pass

        curr_cnf = cnf_out

    print()
    print(f"Done. CSV written: {csv_path}")
    print("Tip: view nicely with:")
    print(f"  column -s, -t {csv_path} | less -S")


if __name__ == "__main__":
    main()
