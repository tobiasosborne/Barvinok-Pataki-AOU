# Meeting Summary: Barvinok--Pataki Theorem and its Generalisation to AOU Spaces

**Date:** 3 March 2026, approx. 10:30--12:00
**Location:** Seminar room (Leibniz Universität Hannover), with Zoom recording
**Duration:** ~82 minutes

## Participants

All participants are members of the Quantum Information groups at the
Institute of Theoretical Physics, Leibniz Universität Hannover.

### In the room

1. **Tobias J. Osborne** (Professor, AG Osborne) — Also logged into Zoom to
   record the session.  Drove the discussion towards generalisation: pushed back
   on basis-dependent arguments ("it's a pain in the ass to choose a basis in
   formalisation"), connected the argument to C\*-algebra state spaces and
   face--commutant correspondences, and led the writing of the AOU space proof
   on the lower blackboard.

2. **Timo Ziegler** (PhD student, AG Osborne; thesis: "Quantum Conic
   Programming") — Main presenter for the first part.  Wrote the classical PSD
   matrix proof on the upper blackboard.  Has detailed notes on the
   Barvinok--Pataki proof and the constructive perturbation argument.

3. **Martin Plavala** (Postdoc, AG Raussendorf; research: general probabilistic
   theories) — Had previously worked out a proof for arbitrary cones
   ("Martin had the proof in general", line ~1406 of the transcript).
   Contributed the general cone perspective to the discussion.

4. **René Schwonnek** (Postdoc / Topical Group Leader for "Quantum Computation
   Concepts", AG Osborne; QuantumFrontiers Cluster of Excellence) — Active in
   the technical debate, particularly the eigenvalue perturbation argument and
   the determinant/characteristic polynomial discussion.

5. **Lennart Binkowski** (PhD, AG Osborne; thesis: "Perfect Tensors from
   Multiple Angles", defended 2025) — Participated in the discussion.

6. **Lukas Hantzko** (PhD student, AG Raussendorf; research: sub-Riemannian
   geometry and measurement-based quantum computation) — Participated in the
   discussion.

7. **Henrik Wilming** (Postdoc, AG Osborne) — Participated in the
   discussion.

### On Zoom

8. **Alexander Stottmeister** (Group Leader, ITP Hannover; research: operator
   algebraic renormalisation) — Asked clarifying questions about the face
   definition (open vs closed condition, extreme rays vs extreme points),
   pointedness of cones, and the infinite-dimensional case.  Flagged the issue
   that the face definition via majorisation looks like an open condition, which
   was resolved (extreme rays, not extreme points).  Raised the question of
   whether the argument needs the commutant to guarantee finite-dimensionality
   of faces.  Had to leave early (daughter waking up).

## Structure of the Discussion

### Part 1: Classical Barvinok--Pataki for PSD Matrices (~0--25 min)

Timo presented a self-contained proof of the rank bound for extreme points of
spectrahedra (feasible sets of semidefinite programs).

**Setup.** Consider the constrained set
$$C = \bigl\{ X \in \mathrm{PSD}_n \;\big|\; \langle A_j, X \rangle = \beta_j,\; j=1,\dots,m \bigr\}.$$
Assume $X \in C$ has $\operatorname{rank}(X) = r$.  Write $X = V B V^\top$ where
$V \in \mathbb{R}^{n \times r}$ and $B \succ 0$ (positive definite, $r \times r$).

**Compression.** Define $\hat{A}_j = V^\top A_j V$ (compressed constraints in
$\mathrm{Sym}_r$) and the linear map
$$\varphi \colon \mathrm{Sym}_r \to \mathbb{R}^m, \qquad
  [\varphi(D)]_j = \langle \hat{A}_j, D \rangle.$$

**Rank--Nullity.** $\dim(\mathrm{Sym}_r) = \frac{r(r+1)}{2}$ and
$\dim(\ker \varphi) = \frac{r(r+1)}{2} - \dim(\operatorname{im} \varphi)
  \geq \frac{r(r+1)}{2} - m$.

**Perturbation.** If $\ker\varphi \neq \{0\}$, pick $D \neq 0$ in
$\ker\varphi$.  Then $B(t) = B + tD$ satisfies all constraints for every
$t \in \mathbb{R}$.  By choosing $t$ so that the smallest eigenvalue of $B(t)$
hits zero (eigenvalue continuity / characteristic polynomial argument), one
obtains a new feasible $B(t_0) \succeq 0$ with
$\operatorname{rank}(B(t_0)) < r$.

**Conclusion.** Iterating, the rank can be reduced until
$\frac{r(r+1)}{2} \leq m$, giving the Barvinok--Pataki bound.  In the complex
case, $\mathrm{Sym}_r$ is replaced by Hermitian matrices of dimension $r^2$,
giving $r^2 \leq m$.

### Part 2: Extended Discussion and Debate (~25--55 min)

Several topics were debated:

- **The perturbation argument in detail.** Significant discussion about *why*
  one can always find a $t$ making $B(t)$ singular.  The key insight: $D$ can
  have negative eigenvalues; you walk $t$ in the direction that pushes the
  smallest eigenvalue of $B + tD$ towards zero.  This is a continuity argument
  (eigenvalues are continuous in $t$), not a "sufficiently small $t$" argument.

- **Feasibility vs optimality.** The optimality case (minimise $\langle C, X
  \rangle$ over $C$) reduces to feasibility with one additional constraint:
  if $D \in \ker\varphi$ but $\langle C, D \rangle \neq 0$, then the current
  point is not optimal.  This is "standard LP/SDP duality" but was discussed
  explicitly for clarity.

- **Basis-free formulation.** Tobias pushed for removing the explicit basis $V$,
  noting that in formalisation (Lean 4) choosing and discarding a basis is
  painful.  The resolution: the $V B V^\top$ decomposition is really saying
  "restrict to the face of the PSD cone containing $X$", which is isomorphic to
  $\mathrm{PSD}_r$.  This face-centric viewpoint generalises directly.

- **Physical applications.** Discussion of what "bite" the theorem has: for
  quantum mechanics with $\leq 3$ constraints (identity + Hamiltonian +
  angular momentum), all extreme points are pure states ($r = 1$).  This means
  pure-state variational methods (e.g.\ Monte Carlo wave function) are still
  valid.  Beyond 3 constraints, mixed states may be needed.

- **One-page PRL.** A half-serious discussion about whether the general
  theorem could fit in a one-page Physical Review Letters paper.

### Part 3: Generalisation to AOU Spaces (~55--82 min)

The group worked together to write out the proof for general Archimedean order
unit spaces on the lower blackboard.  Tobias led the writing.

**Key points of the general proof:**

1. **Setting.** $(V, V^+, e)$ an Archimedean order unit space (AOU space).
   The cone $V^+$ is proper (pointed: $V^+ \cap (-V^+) = \{0\}$) and
   Archimedean (for every $v$, there exists $n$ with $ne \geq v$).

2. **Constrained set.** Given functionals $\alpha_1, \dots, \alpha_m \in V^*$
   and values $b_1, \dots, b_m$:
   $$C = \bigl\{ x \in V^+ \;\big|\; \alpha_j(x) = b_j,\; j = 1,\dots,m \bigr\}.$$

3. **Face of a feasible point.** For $x \in C$, the face
   $F = \operatorname{face}_{V^+}(x)$ is defined as
   $$F = \{ w \in V^+ \mid \exists\, \lambda > 0 :\; x \geq \lambda w \},$$
   i.e.\ all elements of $V^+$ majorised by a positive multiple of $x$.
   Discussion confirmed this is the *smallest face containing $x$ in its
   relative interior*.

4. **The linear map.** Define $\varphi \colon V \to \mathbb{R}^m$ by
   $[\varphi(v)]_j = \alpha_j(v)$.  Restrict to the affine span
   $\operatorname{aff}(F)$ (or equivalently, the linear span $W_F = \operatorname{span}(F)$
   which has dimension $d_F := \dim(F) + 1$ since $F$ is a cone face).

5. **Rank--Nullity on the face.** The restricted map
   $\varphi|_{W_F} \colon W_F \to \mathbb{R}^m$ satisfies:
   $$d_F = \dim(\ker(\varphi|_{W_F})) + \dim(\operatorname{im}(\varphi|_{W_F})) \leq \dim(\ker(\varphi|_{W_F})) + m.$$

6. **Pointed cone → non-trivial kernel → rank reduction.** If
   $\ker(\varphi|_{W_F}) \neq \{0\}$, pick $D \neq 0$ in the kernel.  Because
   the cone is pointed ($V^+ \cap (-V^+) = \{0\}$), the line $\{x + tD : t \in \mathbb{R}\}$
   must eventually leave $V^+$ (in both directions, since otherwise $D$ and
   $-D$ would both be in $V^+$, contradicting pointedness).  Walking along
   this line until hitting a boundary face gives a new feasible point on a
   *strictly smaller face*.

7. **Iteration and bound.** Repeating reduces the face dimension each time.
   The process terminates when $\ker(\varphi|_{W_F}) = \{0\}$, giving:
   $$d_F \leq m \qquad\text{i.e.}\qquad \dim(\operatorname{face}(x)) \leq m - 1.$$

8. **Connection to operator algebras.** For $V = A_{\mathrm{sa}}$ (self-adjoint
   part of a C\*-algebra $A$), the face dimension equals $\dim_{\mathbb{R}}(\pi_\varphi(A)')_{\mathrm{sa}} - 1$
   via the Alfsen--Shultz face--commutant correspondence, recovering
   $\sum n_j^2 \leq m + 1$ (where the $+1$ accounts for the normalisation
   constraint $\varphi(e) = 1$).

### Open Questions and Issues Raised

- **Infinite-dimensional faces.** If the face $F$ is infinite-dimensional, the
  iterative argument becomes an infinite regress.  Alex raised the question of
  whether the commutant is needed to guarantee finiteness.  Resolution: in the
  C\*-algebra case, $m$ constraints force the commutant (and hence the face) to
  be finite-dimensional; in general AOU spaces, finite-dimensionality of the
  face must be assumed or derived from additional structure.

- **Zorn's lemma in infinite dimensions.** The "walk to the boundary" step uses
  supporting hyperplane arguments, which in infinite dimensions invoke Zorn's
  lemma (Hahn--Banach).  In separable spaces, this can be avoided.

- **Face definition subtlety.** The face $\operatorname{face}(x)$ defined via
  majorisation is an *open* condition (strict inequality), but this is
  consistent because extreme rays (not points) are the atoms of cones.

## Action Items

- Tobias to auto-transcribe the recording and produce a LaTeX write-up from
  it (using Claude).
- Timo to share his written proof notes.
- Explore whether the argument can be made fully constructive (no Zorn's lemma)
  in the separable case.
- Consider the "one-page PRL" format for publication.
- Investigate the face dimension gap phenomenon for general AOU spaces.
