#!/usr/bin/env python3
import sys
import csv
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
    # Returns (status, conflicts, cpu_time)
    status = None
    conflicts = None
    cpu_time = None

    for line in text.splitlines():
        s = line.strip()

        # Status line: "s SATISFIABLE" or "s UNSATISFIABLE"
        if s.startswith("s "):
            if "UNSATISFIABLE" in s:
                status = "UNSAT"
            elif "SATISFIABLE" in s:
                status = "SAT"

        # Conflicts line: "c conflicts             : 1301            (....)"
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

        # CPU line: "c CPU time              : 0.008382 s"
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

def main():
    # Usage:
    # python3 versions/asr_v1_active/run_asr_v1.py <input_cnf> <out_prefix> [k_total] [k_step]
    if len(sys.argv) < 3:
        print("Usage: python3 versions/asr_v1_active/run_asr_v1.py <input_cnf> <out_prefix> [k_total] [k_step]")
        print("Example: python3 versions/asr_v1_active/run_asr_v1.py shared/cnf/n120_s2.cnf n120_s2_v1 10 1")
        sys.exit(1)

    input_cnf = Path(sys.argv[1]).expanduser()
    out_prefix = sys.argv[2]
    k_total = int(sys.argv[3]) if len(sys.argv) >= 4 else 10
    k_step  = int(sys.argv[4]) if len(sys.argv) >= 5 else 1

    if k_total <= 0:
        die("k_total must be > 0")
    if k_step <= 0:
        die("k_step must be > 0")
    if k_step > k_total:
        die("k_step cannot be greater than k_total")

    glucose_bin = Path("~/glucose/simp/glucose").expanduser()
    pick_script = Path("versions/asr_v0b/pick_literals_asr_v0b.py")
    apply_units = Path("shared/scripts/apply_units.py")

    if not input_cnf.exists():
        die(f"Input CNF not found: {input_cnf}")
    if not glucose_bin.exists():
        die(f"Glucose binary not found: {glucose_bin}")
    if not pick_script.exists():
        die(f"Pick script not found: {pick_script}")
    if not apply_units.exists():
        die(f"apply_units.py not found: {apply_units}")

    cnf_dir = Path("shared/cnf")
    log_dir = Path("shared/logs")
    res_dir = Path("shared/results")
    cnf_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    res_dir.mkdir(parents=True, exist_ok=True)

    csv_path = res_dir / f"{out_prefix}_asr_v1.csv"

    def solve_and_log(cnf_path: Path, tag: str):
        log_path = log_dir / f"{out_prefix}_{tag}.log"
        rc, out = run_cmd([str(glucose_bin), str(cnf_path)])
        log_path.write_text(out)
        status, conflicts, cpu = parse_glucose_output(out)
        return status, conflicts, cpu, log_path

    def pick_literals(cnf_path: Path, k: int):
        rc, out = run_cmd(["python3", str(pick_script), str(cnf_path), str(k)])
        if rc != 0:
            die(f"pick_literals failed (rc={rc}). Output:\n{out}")

        lits = []
        for tok in out.strip().split():
            try:
                lits.append(int(tok))
            except ValueError:
                pass

        if len(lits) < k:
            die(f"pick_literals returned {len(lits)} literals, expected {k}. Raw output:\n{out}")
        return lits[:k]

    def apply_new_units(curr_cnf: Path, next_cnf: Path, new_lits):
        cmd = ["python3", str(apply_units), str(curr_cnf), str(next_cnf)] + [str(x) for x in new_lits]
        rc, out = run_cmd(cmd)
        if rc != 0:
            die(f"apply_units failed (rc={rc}). Output:\n{out}")
        if not next_cnf.exists():
            die(f"Expected output CNF not created: {next_cnf}")

    # Make a local copy as iter0 CNF under shared/cnf
    curr = cnf_dir / f"{out_prefix}_iter0.cnf"
    curr.write_text(input_cnf.read_text())

    # Write CSV header
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iter", "k_added_total", "new_units", "status", "conflicts", "cpu_time", "cnf_path", "log_path"])

    print("ASR v1 ACTIVE RUN")
    print(f"  input   : {input_cnf}")
    print(f"  prefix  : {out_prefix}")
    print(f"  k_total : {k_total}")
    print(f"  k_step  : {k_step}")
    print()

    # Baseline solve
    status0, conf0, cpu0, log0 = solve_and_log(curr, "iter0")
    with csv_path.open("a", newline="") as f:
        w = csv.writer(f)
        w.writerow([0, 0, "", status0, conf0, cpu0, str(curr), str(log0)])

    print(f"[iter 0] status={status0} conflicts={conf0} cpu={cpu0}  (baseline)")

    k_added = 0
    it = 0

    while k_added < k_total:
        it += 1
        remaining = k_total - k_added
        this_k = k_step if remaining >= k_step else remaining

        # Pick literals on CURRENT CNF
        new_lits = pick_literals(curr, this_k)

        # Apply them to create NEXT CNF
        next_cnf = cnf_dir / f"{out_prefix}_iter{it}.cnf"
        apply_new_units(curr, next_cnf, new_lits)

        k_added += this_k

        # Solve and log
        status, conflicts, cpu, logp = solve_and_log(next_cnf, f"iter{it}")

        # CSV row
        with csv_path.open("a", newline="") as f:
            w = csv.writer(f)
            w.writerow([it, k_added, " ".join(str(x) for x in new_lits), status, conflicts, cpu, str(next_cnf), str(logp)])

        print(f"[iter {it}] +{this_k} units (total={k_added}): {new_lits} | status={status} conflicts={conflicts} cpu={cpu}")

        curr = next_cnf

        # Stop ONLY when solved with 0 conflicts
        if conflicts == 0 and status in ("SAT", "UNSAT"):
            print("Solved with 0 conflicts. Stopping.")
            break

    print()
    print(f"Done. CSV written: {csv_path}")

if __name__ == "__main__":
    main()
