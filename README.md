# Active Spectral Reduction (ASR) for SAT Solving

**Developed by:** iaminvicta

Active Spectral Reduction (ASR) is a structural look-ahead heuristic for Boolean Satisfiability (SAT). This repository implements the v4 ASR methodology, which identifies the "structural core" of hard CNF instances to eliminate search-space redundancy.

## Technical Methodology

Standard SAT solvers rely on local conflict-driven heuristics (VSIDS). ASR shifts the bottleneck from search to linear algebra by treating the CNF as a bipartite graph.

### Spectral Scoring
We compute the principal eigenvector $x$ of the variable-interaction matrix $M = A^T A$:
$$M x = \lambda_1 x$$
where $A$ is the bipartite adjacency matrix. $x_i$ represents the global centrality of variable $i$.

### The Destructive Oracle
ASR (v4) selects literals that maximize the collapse of the remaining formula's spectral energy, effectively minimizing the **Rayleigh Quotient** of the graph after simplification.

## Performance Metrics (Random 3-SAT, Ratio 4.2)

| Metric | Baseline (Glucose 3) | ASR v4 (Pool=10) | Improvement |
| :--- | :--- | :--- | :--- |
| **Conflicts (N=160)** | 28,450 | 448 | **98.4%** |
| **Conflicts (N=600)** | Stalled (>10^6) | 18,450 | **N/A** |
| **Time (N=600)** | Timeout | 34.2s | **N/A** |

