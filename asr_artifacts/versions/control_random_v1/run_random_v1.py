#!/usr/bin/env python3
import sys
import csv
import random
import subprocess
from pathlib import Path

def die(msg: str, code: int = 1):
    print("ERROR:", msg)
    sys.exit(code)

def run_cmd(cmd_list):
    try:
        p = subprocess.run(
            cmd_list,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False
        )
        return p.returncode, p.stdout
    except FileNotFoundError:
        die(f"Command not found: {cmd_list[0]}")

def parse_glucose_output(text: str):
    status = None
    conflicts = None
    cpu_time = None

    for line in text.splitlines():
        s = line.strip()

        if s.startswith("s "):
            if "UNSATISFIABLE" in s:
                status = "UNSAT"
            elif "SATISFIABLE" in s:
                status = "SAT"

        if s.startswith("c conflicts"):
            parts = s.split(":")
            if len(parts) >= 2:
                right = parts[1].strip()
                num = ""
                for ch in right:
                    if ch.isdigit():
                        num += ch
                    else:
                        break
                if num:
                    conflicts = int(num)

        if s.startswith("c CPU time"):
            parts = s.split(":")
            if len(parts) >= 2:
                right = parts[1].strip()
                tok = right.split()
                if tok:
                    try:
                        cpu_time = float(tok[0])
                    except ValueError:
                        pass

    return status, conflicts, cpu_time

def parse_cnf_num_vars(cnf_path: Path):
    with cnf_path.open("r") as f:
        for line in f:
            if line.startswith("p cnf"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    return int(parts[2])
                break
    return None

def main():
    # Usage:
    # python3 versions/control_random_v1/run_random_v1.py <input_cnf> <out_prefix> [k_total] [k_step] [seed]
    if len(sys.argv) < 3:
        print("Usage: python3 versions/control_random_v1/run_random_v1.py <input_cnf> <out_prefix> [k_total] [k_step] [seed]")
        print("Example: python3 versions/control_random_v1/run_random_v1.py shared/cnf/n120_s2.cnf n120_s2_rand 30 1 12345")
        sys.exit(1)

    input_cnf = Path(sys.argv[1]).expanduser()
    out_prefix = sys.argv[2]
    k_total = int(sys.argv[3]) if len(sys.argv) >= 4 else 30
    k_step  = int(sys.argv[4]) if len(sys.argv) >= 5 else 1
    rng_seed = int(sys.argv[5]) if len(sys.argv) >= 6 else 12345

    if k_total <= 0:
        die("k_total must be > 0")
    if k_step <= 0:
        die("k_step must be > 0")
    if k_step > k_total:
        die("k_step cannot be greater than k_total")

    random.seed(rng_seed)

    glucose_bin = Path("~/glucose/simp/glucose").expanduser()
    apply_units = Path("shared/scripts/apply_units.py")

    if not input_cnf.exists():
        die(f"Input CNF not found: {input_cnf}")
    if not glucose_bin.exists():
        die(f"Glucose binary not found: {glucose_bin}")
    if not apply_units.exists():
        die(f"apply_units.py not found: {apply_units}")

    cnf_dir = Path("shared/cnf")
    log_dir = Path("shared/logs")
    res_dir = Path("shared/results")
    cnf_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)

    csv_path = res_dir / f"{out_prefix}_random_v1.csv"

    def solve_and_log(cnf_path: Path, tag: str):
        log_path = log_dir / f"{out_prefix}_{tag}.log"
        rc, out = run_cmd([str(glucose_bin), str(cnf_path)])
        log_path.write_text(out)
        status, conflicts, cpu = parse_glucose_output(out)
        return status, conflicts, cpu, log_path

    def apply_new_units(curr_cnf: Path, next_cnf: Path, new_lits):
        cmd = ["python3", str(apply_units), str(curr_cnf), str(next_cnf)] + [str(x) for x in new_lits]
        rc, out = run_cmd(cmd)
        if rc != 0:
            die(f"apply_units failed (rc={rc}). Output:\n{out}")
        if not next_cnf.exists():
            die(f"Expected output CNF not created: {next_cnf}")

    def pick_random_literals(curr_cnf: Path, want_k: int, used_set):
        # We treat exact literal repeats as forbidden (same as ASR v3 style).
        # Fallback: choose unused variables with positive polarity if needed.
        nvars = parse_cnf_num_vars(curr_cnf)
        if nvars is None:
            die("Could not parse CNF header (missing 'p cnf N M').")

        new_lits = []
        tries = 0

        while len(new_lits) < want_k:
            tries += 1
            if tries > 2000:
                # Fallback: unused vars, + polarity
                for v in range(1, nvars + 1):
                    if len(new_lits) >= want_k:
                        break
                    if v not in used_set and (-v) not in used_set:
                        used_set.add(v)
                        new_lits.append(v)
                return new_lits

            v = random.randint(1, nvars)
            lit = v if random.getrandbits(1) == 0 else -v

            if lit not in used_set:
                used_set.add(lit)
                new_lits.append(lit)

        return new_lits

    # Make iter0 copy in shared/cnf
    curr = cnf_dir / f"{out_prefix}_iter0.cnf"
    curr.write_text(input_cnf.read_text())

    # Write CSV header
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iter", "k_added_total", "new_units", "status", "conflicts", "cpu_time", "cnf_path", "log_path", "rng_seed"])

    print("RANDOM CONTROL v1 (same pipeline, random units)")
    print(f"  input   : {input_cnf}")
    print(f"  prefix  : {out_prefix}")
    print(f"  k_total : {k_total}")
    print(f"  k_step  : {k_step}")
    print(f"  seed    : {rng_seed}")
    print()

    # Baseline solve
    status0, conf0, cpu0, log0 = solve_and_log(curr, "iter0")
    with csv_path.open("a", newline="") as f:
        w = csv.writer(f)
        w.writerow([0, 0, "", status0, conf0, cpu0, str(curr), str(log0), rng_seed])
    print(f"[iter 0] status={status0} conflicts={conf0} cpu={cpu0}  (baseline)")

    used = set()
    k_added = 0
    it = 0

    while k_added < k_total:
        it += 1
        remaining = k_total - k_added
        this_k = k_step if remaining >= k_step else remaining

        new_lits = pick_random_literals(curr, this_k, used)

        next_cnf = cnf_dir / f"{out_prefix}_iter{it}.cnf"
        apply_new_units(curr, next_cnf, new_lits)

        k_added += this_k

        status, conflicts, cpu, logp = solve_and_log(next_cnf, f"iter{it}")

        with csv_path.open("a", newline="") as f:
            w = csv.writer(f)
            w.writerow([it, k_added, " ".join(str(x) for x in new_lits), status, conflicts, cpu, str(next_cnf), str(logp), rng_seed])

        print(f"[iter {it}] +{this_k} units (total={k_added}): {new_lits} | status={status} conflicts={conflicts} cpu={cpu}")

        curr = next_cnf

        if conflicts == 0 and status in ("SAT", "UNSAT"):
            print("Solved with 0 conflicts. Stopping.")
            break

    print()
    print(f"Done. CSV written: {csv_path}")

if __name__ == "__main__":
    main()
