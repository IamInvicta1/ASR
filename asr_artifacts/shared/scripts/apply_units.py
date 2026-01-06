import sys

# Usage:
#   python3 apply_units.py IN.cnf OUT.cnf lit1 lit2 lit3 ...
#
# Example:
#   python3 apply_units.py test_n80.cnf test_n80_units.cnf 1 -2 7

if len(sys.argv) < 4:
    print("usage: python3 apply_units.py IN.cnf OUT.cnf lit1 [lit2 ...]")
    sys.exit(1)

in_cnf = sys.argv[1]
out_cnf = sys.argv[2]
lits = [int(x) for x in sys.argv[3:]]

# Read input CNF and find header
lines = []
num_vars = None
num_clauses = None

with open(in_cnf, "r") as f:
    for line in f:
        s = line.strip()
        if s == "" or s.startswith("c"):
            lines.append(line)
            continue
        if s.startswith("p "):
            parts = s.split()
            # p cnf <vars> <clauses>
            num_vars = int(parts[2])
            num_clauses = int(parts[3])
            lines.append(line)
        else:
            lines.append(line)

if num_vars is None or num_clauses is None:
    print("error: DIMACS header 'p cnf ...' not found")
    sys.exit(1)

# Update clause count in header (add one clause per lit)
new_num_clauses = num_clauses + len(lits)

out_lines = []
header_done = False
for line in lines:
    s = line.strip()
    if (not header_done) and s.startswith("p "):
        out_lines.append(f"p cnf {num_vars} {new_num_clauses}\n")
        header_done = True
    else:
        out_lines.append(line)

# Append unit clauses
for lit in lits:
    out_lines.append(f"{lit} 0\n")

with open(out_cnf, "w") as f:
    f.writelines(out_lines)

print(f"wrote {out_cnf} (added {len(lits)} unit clauses)")
