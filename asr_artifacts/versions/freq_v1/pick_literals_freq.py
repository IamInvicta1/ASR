import sys

# Usage:
#   python3 pick_literals_freq.py IN.cnf K
#
# Output:
#   One line: K literals, space-separated.
#   Heuristic:
#     - pick K variables with highest total occurrences
#     - polarity = sign with majority occurrences (ties -> positive)

def die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)

if len(sys.argv) != 3:
    die("usage: python3 pick_literals_freq.py IN.cnf K")

cnf_path = sys.argv[1]
try:
    K = int(sys.argv[2])
except ValueError:
    die("K must be an integer")

if K <= 0:
    die("K must be > 0")

num_vars = None
pos = {}  # var -> positive count
neg = {}  # var -> negative count

# Parse DIMACS CNF
try:
    with open(cnf_path, "r") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("c"):
                continue

            if s.startswith("p"):
                parts = s.split()
                if len(parts) < 4 or parts[0] != "p" or parts[1] != "cnf":
                    die(f"bad header line: {s}")
                num_vars = int(parts[2])
                # num_clauses = int(parts[3])  # not required for this heuristic
                continue

            # Clause line: integers ending with 0
            # Example: "-12 7 45 0"
            toks = s.split()
            for tok in toks:
                if tok == "0":
                    break
                try:
                    lit = int(tok)
                except ValueError:
                    die(f"non-integer token in clause line: {tok}")

                v = abs(lit)
                if lit > 0:
                    pos[v] = pos.get(v, 0) + 1
                else:
                    neg[v] = neg.get(v, 0) + 1

except FileNotFoundError:
    die(f"file not found: {cnf_path}")
except PermissionError:
    die(f"permission denied reading: {cnf_path}")

if num_vars is None:
    die("DIMACS header 'p cnf <vars> <clauses>' not found")

# Build list of (total_count, var)
counts = []
for v in range(1, num_vars + 1):
    t = pos.get(v, 0) + neg.get(v, 0)
    if t > 0:
        counts.append((t, v))

if not counts:
    # This would mean no clause lines were parsed, which indicates malformed input
    die("no literals found in CNF body (did you generate a valid DIMACS file?)")

# Sort by total occurrences desc, tie by var asc for determinism
counts.sort(key=lambda x: (-x[0], x[1]))

picked = counts[:K]

lits = []
for total, v in picked:
    p = pos.get(v, 0)
    n = neg.get(v, 0)
    lit = v if p >= n else -v  # tie -> positive
    lits.append(str(lit))

print(" ".join(lits))

