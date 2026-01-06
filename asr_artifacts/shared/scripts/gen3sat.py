import random
import sys

# Random 3-SAT generator in DIMACS CNF
# Usage: python3 gen3sat.py N M SEED OUTFILE

if len(sys.argv) != 5:
    print("usage: python3 gen3sat.py N M SEED OUTFILE")
    sys.exit(1)

N = int(sys.argv[1])
M = int(sys.argv[2])
SEED = int(sys.argv[3])
OUTFILE = sys.argv[4]

random.seed(SEED)

with open(OUTFILE, "w") as f:
    f.write(f"p cnf {N} {M}\n")
    for _ in range(M):
        vars3 = random.sample(range(1, N + 1), 3)
        lits = []
        for v in vars3:
            if random.getrandbits(1) == 0:
                lits.append(str(v))
            else:
                lits.append(str(-v))
        f.write(" ".join(lits) + " 0\n")
