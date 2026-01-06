import sys
import math

# ASR v0b (Option B): incidence-graph spectral scoring (variables <-> clauses)
#
# Usage:
#   python3 pick_literals_asr_v0b.py IN.cnf K
#   python3 pick_literals_asr_v0b.py IN.cnf K ITERS
#
# Output:
#   One line: K literals, space-separated.
#
# Method (single pass):
#   - Parse DIMACS CNF
#   - Build bipartite incidence structure:
#       clauses: list of literals
#       var_to_clauses: for each var, which clauses contain it (ignoring sign)
#   - Power-iterate on bipartite graph:
#       clause_score[j] = sum(var_score[v] for v in vars(clause j))
#       var_score[v]    = sum(clause_score[j] for j in clauses_containing_v)
#     Normalize each step.
#   - Pick top-K vars by var_score
#   - Polarity: majority polarity in original CNF (ties -> positive)

def die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)

def l2_normalize(vec):
    s = 0.0
    for x in vec:
        s += x * x
    if s <= 0.0:
        return vec
    inv = 1.0 / math.sqrt(s)
    return [x * inv for x in vec]

def parse_dimacs(path):
    num_vars = None
    clauses = []
    pos = {}
    neg = {}

    try:
        with open(path, "r") as f:
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
                    continue

                # clause: ints ending with 0
                toks = s.split()
                lits = []
                for tok in toks:
                    if tok == "0":
                        break
                    lit = int(tok)
                    v = abs(lit)
                    if lit > 0:
                        pos[v] = pos.get(v, 0) + 1
                    else:
                        neg[v] = neg.get(v, 0) + 1
                    lits.append(lit)

                # ignore empty clauses (shouldn't happen in your generator output)
                if lits:
                    clauses.append(lits)

    except FileNotFoundError:
        die(f"file not found: {path}")
    except PermissionError:
        die(f"permission denied reading: {path}")

    if num_vars is None:
        die("DIMACS header 'p cnf <vars> <clauses>' not found")
    if not clauses:
        die("no clauses parsed from CNF body")

    return num_vars, clauses, pos, neg

def build_incidence(num_vars, clauses):
    # var_to_clauses[v] = list of clause indices containing v (v is 1..num_vars)
    var_to_clauses = [[] for _ in range(num_vars + 1)]
    clause_vars = []  # per clause: unique vars (ignore sign)
    for j, lits in enumerate(clauses):
        seen = set()
        vars_j = []
        for lit in lits:
            v = abs(lit)
            if v not in seen:
                seen.add(v)
                vars_j.append(v)
                var_to_clauses[v].append(j)
        clause_vars.append(vars_j)
    return var_to_clauses, clause_vars

def spectral_scores(num_vars, var_to_clauses, clause_vars, iters):
    m = len(clause_vars)

    # start with uniform var scores on vars that actually appear
    var_score = [0.0] * (num_vars + 1)
    active = 0
    for v in range(1, num_vars + 1):
        if var_to_clauses[v]:
            var_score[v] = 1.0
            active += 1

    if active == 0:
        die("no active variables found")

    # normalize initial
    vs = [var_score[v] for v in range(1, num_vars + 1)]
    vs = l2_normalize(vs)
    for v in range(1, num_vars + 1):
        var_score[v] = vs[v - 1]

    clause_score = [0.0] * m

    for _ in range(iters):
        # clause scores from vars
        for j in range(m):
            s = 0.0
            for v in clause_vars[j]:
                s += var_score[v]
            clause_score[j] = s

        # normalize clause scores
        clause_score = l2_normalize(clause_score)

        # var scores from clauses
        new_var = [0.0] * (num_vars + 1)
        for v in range(1, num_vars + 1):
            s = 0.0
            for j in var_to_clauses[v]:
                s += clause_score[j]
            new_var[v] = s

        # normalize var scores (excluding index 0)
        vs = [new_var[v] for v in range(1, num_vars + 1)]
        vs = l2_normalize(vs)
        for v in range(1, num_vars + 1):
            var_score[v] = vs[v - 1]

    return var_score

def pick_top_k_literals(num_vars, var_score, pos, neg, k):
    scored = []
    for v in range(1, num_vars + 1):
        s = var_score[v]
        # skip vars that never appeared (score likely 0)
        if s > 0.0:
            scored.append((s, v))

    if not scored:
        die("no variables had positive score (unexpected)")

    scored.sort(key=lambda x: (-x[0], x[1]))
    picked = scored[:k]

    lits = []
    for s, v in picked:
        p = pos.get(v, 0)
        n = neg.get(v, 0)
        lit = v if p >= n else -v
        lits.append(str(lit))
    return " ".join(lits)

def main():
    if len(sys.argv) not in (3, 4):
        die("usage: python3 pick_literals_asr_v0b.py IN.cnf K [ITERS]")

    cnf_path = sys.argv[1]
    try:
        k = int(sys.argv[2])
    except ValueError:
        die("K must be an integer")
    if k <= 0:
        die("K must be > 0")

    iters = 20
    if len(sys.argv) == 4:
        try:
            iters = int(sys.argv[3])
        except ValueError:
            die("ITERS must be an integer")
        if iters <= 0:
            die("ITERS must be > 0")

    num_vars, clauses, pos, neg = parse_dimacs(cnf_path)
    var_to_clauses, clause_vars = build_incidence(num_vars, clauses)
    var_score = spectral_scores(num_vars, var_to_clauses, clause_vars, iters)
    out = pick_top_k_literals(num_vars, var_score, pos, neg, k)
    print(out)

if __name__ == "__main__":
    main()

