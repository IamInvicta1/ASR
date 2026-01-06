#!/usr/bin/env python3
import os
import sys
import csv
import subprocess

# ----------------------------
# DIMACS CNF helpers
# ----------------------------
def read_dimacs(path):
    comments = []
    nvars = 0
    clauses = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("c"):
                comments.append(line)
                continue
            if line.startswith("p"):
                parts = line.split()
                if len(parts) >= 4 and parts[1] == "cnf":
                    nvars = int(parts[2])
                continue
            parts = line.split()
            lits = [int(x) for x in parts]
            if not lits:
                continue
            if lits[-1] == 0:
                lits = lits[:-1]
            clauses.append(lits)
    if nvars <= 0:
        raise ValueError("Missing/invalid DIMACS header (p cnf ...).")
    return nvars, clauses, comments

def write_dimacs(path, nvars, clauses, comments):
    with open(path, "w") as f:
        for c in comments:
            f.write(c + "\n")
        f.write("p cnf {} {}\n".format(nvars, len(clauses)))
        for cl in clauses:
            if len(cl) == 0:
                f.write("0\n")
            else:
                f.write(" ".join(str(x) for x in cl) + " 0\n")

def add_unit_clauses(src_cnf, dst_cnf, new_units):
    nvars, clauses, comments = read_dimacs(src_cnf)
    for lit in new_units:
        clauses.append([int(lit)])
    write_dimacs(dst_cnf, nvars, clauses, comments)

# ----------------------------
# Glucose runner + log parser
# ----------------------------
def run_glucose(glucose_bin, cnf_path, log_path):
    with open(log_path, "w") as f:
        subprocess.run([glucose_bin, cnf_path], stdout=f, stderr=subprocess.STDOUT, check=False)

def parse_glucose_log(log_path):
    status = "UNKNOWN"
    conflicts = 10**18
    cpu_time = 0.0
    with open(log_path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("s "):
                if "SATISFIABLE" in s and "UNSAT" not in s:
                    status = "SAT"
                elif "UNSATISFIABLE" in s:
                    status = "UNSAT"
            if s.startswith("c conflicts"):
                parts = s.replace(":", " ").split()
                for i in range(len(parts)):
                    if parts[i] == "conflicts" and i + 1 < len(parts):
                        try:
                            conflicts = int(parts[i+1])
                        except:
                            pass
            if s.startswith("c CPU time"):
                parts = s.replace(":", " ").split()
                try:
                    cpu_time = float(parts[-2]) if parts[-1] == "s" else float(parts[-1])
                except:
                    pass
    return status, conflicts, cpu_time

# ----------------------------
# V4 proxy WITHOUT numpy
# power-iteration on bipartite walk
# ----------------------------
def build_clause_vars(clauses, nvars):
    clause_vars = []
    for cl in clauses:
        vs = []
        for lit in cl:
            v = abs(lit)
            if 1 <= v <= nvars:
                vs.append(v)
        clause_vars.append(vs)
    return clause_vars

def power_iter_var_scores(clause_vars, nvars, iters=25):
    x = [0.0] * (nvars + 1)
    for v in range(1, nvars + 1):
        x[v] = 1.0

    for _ in range(iters):
        y = [0.0] * len(clause_vars)
        for j, vs in enumerate(clause_vars):
            s = 0.0
            for v in vs:
                s += x[v]
            y[j] = s

        x2 = [0.0] * (nvars + 1)
        for j, vs in enumerate(clause_vars):
            val = y[j]
            for v in vs:
                x2[v] += val

        norm = 0.0
        for v in range(1, nvars + 1):
            norm += x2[v] * x2[v]
        norm = norm ** 0.5
        if norm == 0.0:
            break
        for v in range(1, nvars + 1):
            x[v] = x2[v] / norm

    return x

def literal_scores_and_polarity(clauses, nvars, iters=25):
    clause_vars = build_clause_vars(clauses, nvars)
    x = power_iter_var_scores(clause_vars, nvars, iters=iters)

    clause_mass = [0.0] * len(clauses)
    for j, cl in enumerate(clauses):
        m = 0.0
        for lit in cl:
            m += x[abs(lit)]
        clause_mass[j] = m

    pos = [0.0] * (nvars + 1)
    neg = [0.0] * (nvars + 1)
    for j, cl in enumerate(clauses):
        w = clause_mass[j]
        for lit in cl:
            v = abs(lit)
            if lit > 0:
                pos[v] += w
            else:
                neg[v] += w

    lit_score = {}
    for v in range(1, nvars + 1):
        strength = x[v] * (pos[v] + neg[v] + 1e-12)
        pol = 1 if pos[v] >= neg[v] else -1
        lit_score[v] = (strength, pol)
    return lit_score

def pick_units_rayleigh(cnf_path, k, iters=25, avoid=set()):
    nvars, clauses, _ = read_dimacs(cnf_path)
    scores = literal_scores_and_polarity(clauses, nvars, iters=iters)
    order = list(range(1, nvars + 1))
    order.sort(key=lambda v: scores[v][0], reverse=True)

    chosen = []
    for v in order:
        pol = scores[v][1]
        lit = pol * v
        if lit in avoid or -lit in avoid:
            continue
        chosen.append(lit)
        if len(chosen) >= k:
            break
    return chosen

# ----------------------------
# Main: FORCE WALK until conflicts==0 or max_iters
# ----------------------------
def ensure_dir(p):
    if p and not os.path.exists(p):
        os.makedirs(p, exist_ok=True)

def main():
    if len(sys.argv) != 6:
        print("Usage: python3 run_asr_v4_forcewalk.py <input_cnf> <tag> <max_iters> <k_per_iter> <glucose_bin>")
        sys.exit(1)

    input_cnf = sys.argv[1]
    tag = sys.argv[2]
    max_iters = int(sys.argv[3])
    k_per_iter = int(sys.argv[4])
    glucose_bin = sys.argv[5]

    ensure_dir("shared/cnf")
    ensure_dir("shared/logs")
    ensure_dir("shared/results")

    csv_path = "shared/results/{}_v4_forcewalk.csv".format(tag)

    curr = input_cnf
    avoid = set()
    k_total = 0

    with open(csv_path, "w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["iter","k_added_total","new_units","status","conflicts","cpu_time","cnf_path","log_path"])

        # iter 0: solve baseline and log it, but DO NOT STOP
        log0 = "shared/logs/{}_v4_fw_iter0.log".format(tag)
        run_glucose(glucose_bin, curr, log0)
        status0, conf0, cpu0 = parse_glucose_log(log0)
        print("[iter 0] status={} conflicts={} cpu={}".format(status0, conf0, cpu0))
        w.writerow([0,0,"",status0,conf0,cpu0,curr,log0])

        if conf0 == 0 and status0 in ("SAT","UNSAT"):
            print("Already 0 conflicts at baseline. Done.")
            print("Done. CSV written:", csv_path)
            return

        for it in range(1, max_iters + 1):
            new_units = pick_units_rayleigh(curr, k_per_iter, iters=25, avoid=avoid)
            if len(new_units) < k_per_iter:
                print("ERROR: not enough unique units found; try smaller k_per_iter.")
                break

            for lit in new_units:
                avoid.add(lit)
                avoid.add(-lit)

            k_total += len(new_units)

            cnf_out = "shared/cnf/{}_v4_fw_iter{}.cnf".format(tag, it)
            log_out = "shared/logs/{}_v4_fw_iter{}.log".format(tag, it)

            add_unit_clauses(curr, cnf_out, new_units)
            run_glucose(glucose_bin, cnf_out, log_out)
            status, conflicts, cpu = parse_glucose_log(log_out)

            print("[iter {}] +{} units (total={}): {} | status={} conflicts={} cpu={}".format(
                it, len(new_units), k_total, new_units, status, conflicts, cpu
            ))

            w.writerow([it,k_total," ".join(str(x) for x in new_units),status,conflicts,cpu,cnf_out,log_out])

            if conflicts == 0 and status in ("SAT","UNSAT"):
                print("Solved with 0 conflicts. Stopping.")
                break

            curr = cnf_out

    print("Done. CSV written:", csv_path)

if __name__ == "__main__":
    main()
