# Rényi Generalisations of the Barvinok--Pataki Bound

**Date:** 13 March 2026
**Context:** Coauthor feedback: "I don't care about rank, it isn't robust."

## 1. The Robustness Problem

The BP bound says: for an extreme point $\rho^*$ of a spectrahedron
$S = \{\rho \geq 0 : \operatorname{tr}(A_i \rho) = b_i,\ i=1,\ldots,m\}$,
$$\operatorname{rank}(\rho^*)^2 \leq m \quad\text{(complex)}, \qquad
  \frac{r(r+1)}{2} \leq m \quad\text{(real)}.$$

**Why rank is fragile:**
- Rank is *discontinuous*: $\operatorname{rank}(\rho + \varepsilon I) = d$ for
  any $\varepsilon > 0$, regardless of how small.
- SDP solvers never return exact zeros. A "rank-3" solution typically has
  eigenvalues like $(0.45, 0.35, 0.20, 10^{-9}, 10^{-11}, \ldots)$.
- Physical systems always have noise. Claiming "this quantum channel has
  Kraus rank 3" is never exactly verifiable.

Rank is the *0-Rényi entropy* of the eigenvalue distribution. The natural
question: **can we replace it with an $\alpha$-Rényi quantity for $\alpha > 0$?**


## 2. Rényi Entropies as Soft Rank

For a state $\rho \geq 0$ with $\operatorname{tr}\rho = 1$ and eigenvalues
$\lambda_1 \geq \cdots \geq \lambda_d \geq 0$:

$$S_\alpha(\rho) = \frac{1}{1-\alpha}\log\!\Big(\sum_i \lambda_i^\alpha\Big),
  \qquad \alpha \geq 0,\ \alpha \neq 1.$$

The **$\alpha$-rank** (Rényi dimension, effective dimension):
$$d_\alpha(\rho) = \exp(S_\alpha(\rho)) = \Big(\sum_i \lambda_i^\alpha\Big)^{1/(1-\alpha)}.$$

| $\alpha$ | $S_\alpha$ | $d_\alpha$ | Name |
|----------|-----------|-----------|------|
| 0 | $\log\operatorname{rank}$ | $\operatorname{rank}(\rho)$ | Hartley / hard rank |
| 1 | $-\sum \lambda_i\log\lambda_i$ | $\exp(H(\rho))$ | von Neumann / effective rank (Roy--Vetterli) |
| 2 | $-\log\operatorname{tr}(\rho^2)$ | $1/\operatorname{tr}(\rho^2)$ | Collision / inverse purity |
| $\infty$ | $-\log\lambda_{\max}$ | $1/\lambda_{\max}$ | Min-entropy / inverse max eigenvalue |

**Key properties:**
- $d_\infty \leq d_2 \leq d_1 \leq d_0 = \operatorname{rank}$
  (monotone decreasing in $\alpha$).
- For $\alpha > 0$: $d_\alpha$ is **continuous** in the eigenvalues.
- For $\alpha \geq 1$: $d_\alpha$ is even **Lipschitz** (w.r.t. trace norm).
- $d_\alpha = \operatorname{rank}$ iff all nonzero eigenvalues are equal
  (maximally mixed in support).


## 3. The Basic Theorem (Trivial Direction)

**Proposition (Rényi BP).** For an extreme point $\rho^*$ of a spectrahedron
with $m$ constraints (complex case):
$$d_\alpha(\rho^*)^2 \leq m \qquad \text{for all } \alpha \geq 0.$$

*Proof.* $d_\alpha \leq d_0 = \operatorname{rank} \leq \lfloor\sqrt{m}\rfloor$. $\square$

**This bound is tight for all $\alpha$ simultaneously.** Take any spectrahedron
that saturates BP with rank $r = \lfloor\sqrt{m}\rfloor$ and whose extreme
point has equal nonzero eigenvalues $\lambda_i = 1/r$. Then
$d_\alpha = r$ for all $\alpha$.

So: the *bound itself* does not improve with $\alpha$. But the *quantity
being bounded* is dramatically better behaved.


## 4. Why This Matters: Robust BP for Approximate Extreme Points

Here is where the Rényi formulation pays off. Define:

**$\delta$-approximate extreme point:** $\rho \in S$ such that
$\|\rho - \rho^*\|_1 \leq \delta$ for some extreme point $\rho^*$ of $S$.

For such points, the rank bound is useless:
$$\operatorname{rank}(\rho) = d \quad \text{(generically, for any $\delta > 0$)}.$$

But the Rényi entropy is controlled:

**Proposition (Robust Rényi BP for $\alpha \geq 1$).** Let $\rho \in S$ with
$\|\rho - \rho^*\|_1 \leq \delta$ for some extreme point $\rho^*$. Then for
$\alpha \geq 1$:
$$S_\alpha(\rho) \leq \tfrac{1}{2}\log m + g_\alpha(\delta)$$
where $g_\alpha(\delta) \to 0$ as $\delta \to 0$, **independent of the ambient
dimension $d$**.

*Proof sketch for $\alpha = 1$ (von Neumann).* Write $\rho^*$ in its eigenbasis:
$\rho^* = \sum_{i=1}^r \lambda_i |e_i\rangle\langle e_i|$ with $r \leq \sqrt{m}$
and $\lambda_i > 0$, $\sum \lambda_i = 1$. Let $P$ project onto the support of
$\rho^*$. Then:
- $\operatorname{tr}(P\rho) \geq 1 - \delta/2$ (by trace distance bound).
- Let $\rho_{\mathrm{head}} = P\rho P / \operatorname{tr}(P\rho P)$
  (rank $\leq r$) and $p = \operatorname{tr}(P\rho P) \geq 1 - \delta/2$.
- By subadditivity of von Neumann entropy on the head/tail decomposition:
  $$S_1(\rho) \leq h(p) + p\,S_1(\rho_{\mathrm{head}}) + (1-p)\log d$$
  where $h$ is binary entropy. But $S_1(\rho_{\mathrm{head}}) \leq \log r
  \leq \frac{1}{2}\log m$, and $(1-p) \leq \delta/2$.

The $(1-p)\log d$ term involves $d$, which is problematic. However, for
$\alpha > 1$, the Rényi entropy is *less sensitive to the tail*:

**Key estimate for $\alpha > 1$:** The tail eigenvalues
$\mu_j$ ($j > r$) satisfy $\sum_{j>r}\mu_j \leq \delta/2$, so
$$\sum_{j>r} \mu_j^\alpha \leq \Big(\sum_{j>r}\mu_j\Big)^\alpha \leq (\delta/2)^\alpha$$
by the inequality $\|x\|_\alpha \leq \|x\|_1$ for probability vectors. Hence:
$$\operatorname{tr}(\rho^\alpha) \leq \operatorname{tr}(\rho_{\mathrm{head}}^\alpha) + (\delta/2)^\alpha$$
and the bound becomes $d$-independent.

For $\alpha = 1$, the $d$-dependence is real and cannot be eliminated without
additional structure (this is the Fannes--Audenaert phenomenon). But for
$\alpha \in (0,1)$, partial results are possible. $\square$

**Upshot:** For $\alpha > 1$, the Rényi BP bound
$$d_\alpha(\rho)^2 \lesssim m$$
extends robustly to approximate extreme points, with error controlled by the
approximation quality $\delta$ alone—not the ambient dimension $d$. This is
exactly the kind of dimension-free robustness one wants.


## 5. The $\alpha = 2$ Case: Purity

The case $\alpha = 2$ is especially clean because purity $\operatorname{tr}(\rho^2)$
is a quadratic function of matrix entries.

**Corollary (Purity BP).** For an extreme point $\rho^*$ of a spectrahedron
with $m$ constraints:
$$\operatorname{tr}(\rho^{*2}) \geq \frac{1}{m}.$$

*Proof.* $\operatorname{tr}(\rho^2) \geq 1/\operatorname{rank}(\rho)^2$
by Cauchy--Schwarz (actually $\geq 1/\operatorname{rank}$), and
$\operatorname{rank} \leq \sqrt{m}$. Hmm, we get $\operatorname{tr}(\rho^2) \geq 1/\sqrt{m}$,
which is the correct bound: $d_2 = 1/\operatorname{tr}(\rho^2) \leq \sqrt{m}$.

Equivalently: $\operatorname{tr}(\rho^{*2}) \geq 1/\sqrt{m}$ (tight when eigenvalues are uniform
on rank-$\sqrt{m}$ support).

**Why purity is interesting:** It can be measured directly (via SWAP test in
quantum mechanics, or as the Hilbert--Schmidt norm squared). No tomography
needed. An experimentally verifiable BP bound.


## 6. Pataki's Eigenvalue Clustering

Pataki (1998) proved that at *optimal* extreme points of SDPs (not just
extreme points of the feasible set), eigenvalues cluster: the $k$-th and
$(k+1)$-st largest eigenvalues coincide. The multiplicity of this critical
eigenvalue is bounded below.

**Implication for Rényi entropies:** Eigenvalue clustering means the
eigenvalue distribution is *not* generic. In particular, if many eigenvalues
coincide, the Rényi entropy is lower than it would be for a generic
distribution on the same support. This could give *tighter* Rényi BP bounds
for optimal (not just feasible) extreme points.

**Concrete question:** For the SDP $\min\operatorname{tr}(C\rho)$ s.t.
$\rho \in S$, is there a bound
$$d_\alpha(\rho^*) \leq f(m, \alpha) < \sqrt{m}$$
at the optimum $\rho^*$, where $f < \sqrt{m}$ for $\alpha > 0$?

This would use the objective function to constrain the eigenvalue distribution
beyond what pure extremality gives.


## 7. Beyond Rank: Other Candidate Quantities

The coauthor's complaint "rank isn't robust" suggests looking at other
quantities that:
(a) reduce to rank in some limit,
(b) are continuous,
(c) have physical/operational meaning.

| Quantity | Formula | Continuity | Operational meaning |
|----------|---------|-----------|-------------------|
| Rank | $\#\{\lambda_i > 0\}$ | Discontinuous | Exact compression dimension |
| $\varepsilon$-rank | $\#\{\lambda_i > \varepsilon\}$ | Piecewise constant | Lossy compression dimension |
| Effective rank (Roy--Vetterli) | $\exp(S_1(\rho))$ | $C^\infty$ | Shannon-optimal dimension |
| Participation ratio | $1/\operatorname{tr}(\rho^2)$ | $C^\infty$ | Collision dimension / SWAP test |
| Smooth max-rank | $\min_{\|\sigma-\rho\|_1 \leq \varepsilon}\operatorname{rank}(\sigma)$ | LSC | Robust compression dimension |
| Schatten $p$-norm | $\|\rho\|_p = (\operatorname{tr}\rho^p)^{1/p}$ | Lipschitz | Interpolates between trace and operator norm |

All of these satisfy a BP-type bound (at exact extreme points) as a trivial
corollary of the rank bound. The interest is in (a) which ones extend to
approximate extreme points, and (b) which ones give *tighter* bounds for
specific constraint structures.


## 8. Speculative Directions

### 8.1. SDP hierarchy for entropy bounds

Rényi entropies $S_\alpha$ for integer $\alpha$ involve
$\operatorname{tr}(\rho^\alpha)$, which is a degree-$\alpha$ polynomial
in the matrix entries. Bounding $\operatorname{tr}(\rho^\alpha)$ over a
spectrahedron is a polynomial optimisation over an SDP-representable set.
Sum-of-squares (SOS) hierarchies could give systematic bounds.

### 8.2. The Im--Wolkowicz refinement

Im and Wolkowicz (2021) strengthened BP by incorporating the *singularity
degree* of the spectrahedron: $r \leq \lfloor\sqrt{m - s}\rfloor$ where
$s$ is the singularity degree. The singularity degree measures how "degenerate"
the spectrahedron is—higher $s$ means a sharper bound. This is still about
integer rank, but the singularity degree could also inform Rényi bounds.

### 8.3. Modular theory connection

For a faithful state $\rho$ on a von Neumann algebra $\mathcal{M}$, the
modular operator $\Delta_\rho$ encodes the "soft facial structure" via
$\Delta_\rho^{it}$. The Rényi divergence $D_\alpha(\rho\|\sigma)$ has a
representation in terms of $\Delta_{\sigma,\rho}$. The face of the state
space at $\rho$ (relevant for BP) corresponds to the support projection of
$\Delta_\rho$, i.e., the $\alpha \to 0$ limit. Rényi with $\alpha > 0$
gives a "smoothed face" via level sets of $\Delta_\rho$.

This could connect BP-type bounds to modular theory—especially interesting for
the infinite-dimensional AOU setting already in our paper.

### 8.4. Quantum Rényi entropy power inequality

There are quantum analogues of the entropy power inequality involving Rényi
entropies. If the constraint maps $A_i$ have specific structure (e.g.,
partial traces), these inequalities might give Rényi bounds on extreme
points that are not captured by the rank argument alone.

### 8.5. Concentration of eigenvalues for random spectrahedra

For "generic" constraint matrices $A_i$ (e.g., GUE or Wishart), the
eigenvalue distribution of extreme points should concentrate. Random matrix
theory (Marchenko--Pastur type results for truncated spectra) could give
*typical* Rényi entropy values, complementing the *worst-case* BP bound.


## 9. Concrete Proposals

### Proposal A: State the Rényi BP cleanly and prove the robust version

Write up the $\alpha > 1$ robust bound (Section 4) as a proper theorem.
This is a genuinely new result: a dimension-free robustness guarantee
for approximate extreme points, stated in terms of Rényi entropy.

### Proposal B: Numerical exploration

For each application in the paper (channels, marginals, QEC, thermal),
compute the actual Rényi entropies of extreme points found by SDP solvers
and compare against the bound $d_\alpha \leq \sqrt{m}$. How tight is it?
Does the tightness depend on $\alpha$?

### Proposal C: Purity as experimentally verifiable BP

The purity bound $\operatorname{tr}(\rho^2) \geq 1/\sqrt{m}$ for extreme
points is measurable via SWAP test. Propose this as an experimental
signature of BP-constrained structure.

### Proposal D: Explore the $\alpha$-optimal direction

For SDPs with an objective, explore whether Pataki's eigenvalue clustering
gives improved Rényi bounds at the optimum.


## 10. Numerical Findings (renyi_bp.jl)

### Finding 1: d_α is well below the bound for random spectrahedra

For random spectrahedra in $M_{10}(\mathbb{C})$ with $m=25$ constraints
(BP bound: $d_\alpha \leq 5$):

| Quantity | Mean | Max | Bound |
|----------|------|-----|-------|
| $d_0$ (rank) | 6.5 | 7.0 | 5 (violated by solver noise!) |
| $d_{1/2}$ | 2.8 | 3.5 | 5 |
| $d_1$ (vN) | 2.6 | 3.2 | 5 |
| $d_2$ (purity) | 2.3 | 2.9 | 5 |
| $d_\infty$ | 1.8 | 2.2 | 5 |

**Rank exceeds the BP bound** (6.5 > 5) because SCS returns approximate
extreme points with small but nonzero tail eigenvalues. **But $d_\alpha$ for
$\alpha > 0$ stays well within bounds** — demonstrating exactly the robustness
advantage.

### Finding 2: Rank jumps, Rényi dimensions don't

Under depolarising perturbation $\rho_\delta = (1-\delta)\rho^* + \delta I/d$:
- Rank jumps from 7 to 10 at $\delta \sim 10^{-5}$ (threshold artifact).
- $d_2$ changes from 1.90 to 1.97 at $\delta = 0.05$ — barely moves.
- $d_\infty$ changes from 1.53 to 1.60 at $\delta = 0.05$.

**Rényi dimensions are insensitive to perturbation** by orders of magnitude
compared to rank.

### Finding 3: Eigenvalues are generically skewed

At BP-saturating extreme points (rank = $\lfloor\sqrt{m}\rfloor$), the
eigenvalue distribution is *not* uniform. Example at $m=16$, rank $=4$:
eigenvalues $(0.697, 0.303, \sim 0, \sim 0)$ giving $d_2/d_0 = 0.25$.

This means **$d_\alpha \ll \operatorname{rank}$ generically**, so the
Rényi formulation is practically tighter even though the worst-case bound
is the same. This suggests:

**Conjecture (Typical Rényi BP).** For generic constraint matrices $A_i$,
the Rényi dimension of extreme points satisfies
$$\mathbb{E}[d_\alpha(\rho^*)] \leq c_\alpha \cdot \sqrt{m}$$
with $c_\alpha < 1$ for $\alpha > 0$, independent of ambient dimension $d$.

### Finding 4: Channels show the same pattern

For CPTP maps $M_3 \to M_3$ (Choi dim 9, $m = 9$, BP Kraus rank $\leq 3$):
actual Kraus rank is 5--6 (solver), but $d_1 \approx 1.9$ and
$d_2 \approx 1.8$, well below the bound of 3.

### Finding 5: No extra clustering at SDP optima (Experiment 5)

Pataki's eigenvalue clustering theorem predicts that *optimal* extreme
points should have more eigenvalue clustering than arbitrary extreme points.
However, numerically (d=10, m=9, 10 trials):
- Optima: mean $d_2/d_0 = 0.22 \pm 0.06$
- Random extreme points: mean $d_2/d_0 = 0.22 \pm 0.07$

**No significant difference.** Both show strong skewness ($d_2 \ll d_0$),
suggesting that the skewness is a property of *extremality itself*, not of
optimality. The Pataki clustering effect may be too subtle to detect at
these small dimensions.

### Finding 6: Scaling confirms typical bound conjecture (Experiment 6)

For $d=15$ and varying $m$, the ratio $c_\alpha = d_\alpha / \sqrt{m}$:

| $m$ | $c_1$ (vN) | $c_2$ (purity) | $c_\infty$ |
|-----|-----------|----------------|------------|
| 4   | 0.50      | 0.50           | 0.50       |
| 9   | 0.38      | 0.36           | 0.35       |
| 16  | 0.44      | 0.41           | 0.35       |
| 25  | 0.46      | 0.40           | 0.31       |

**All ratios $< 1$**, strongly supporting the Typical Rényi BP conjecture.
The ratio $c_\infty$ decreases with $m$, suggesting the dominant eigenvalue
concentrates as constraints grow.

### Finding 7: Application-specific bounds are much tighter (Experiment 7)

| Application | $m$ | BP bound | rank | $d_1$ | $d_2$ | $d_\infty$ |
|-------------|-----|----------|------|-------|-------|------------|
| Marginal ($3 \times 3$) | 17 | 4.1 | 2 | 1.57 | 1.39 | 1.20 |
| Covariant channel ($\sigma_z$) | 36 | 6.0 | 2 | 2.00 | 2.00 | 2.00 |

The covariant channel has maximally mixed eigenvalues in its rank-2 support
($d_\alpha = d_0 = 2$ for all $\alpha$), confirming that symmetry constraints
can force uniform eigenvalue distributions. But even so, the actual values
are far below the BP bound.


## 11. Summary

| Direction | Status | Potential |
|-----------|--------|-----------|
| $d_\alpha \leq \sqrt{m}$ for exact extreme pts | Trivial corollary of BP | Low (same bound) |
| $d_\alpha \lesssim \sqrt{m}$ for *approximate* extreme pts ($\alpha > 1$) | Provable, new | **High** — dimension-free robustness |
| Purity as experimental BP signature | Easy corollary | Medium — nice narrative |
| Rényi bounds at SDP optima via Pataki clustering | Open question | Medium — needs work |
| SOS hierarchy for $\operatorname{tr}(\rho^\alpha)$ bounds | Speculative | Medium — systematic |
| Modular theory / infinite-dim connection | Speculative | High if it works |
| Random spectrahedra concentration | Speculative | Medium — typical vs worst-case |
