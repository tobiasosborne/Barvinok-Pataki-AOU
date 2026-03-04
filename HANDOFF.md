# Handoff Report: Barvinok-Pataki for Order Unit Spaces

## Status

**Paper at 17 pages with conic inequality extensions, spectrahedra, and Cuntz algebra example.**

- LaTeX write-up: `discussion/barvinok_pataki_aou_proof.tex` (17 pages, compiles cleanly)
- Numerical scripts: `numerics/spectrahedra_bp.jl`, `numerics/cuntz_sdp.jl`
- Email synthesis: `discussion/email_extensions.md`
- 7 bibliography entries, 12 numbered sections + acknowledgements

---

## Work Completed (Session 4)

### Cuntz Algebra O_2 Example (Section 12)

New infinite-dimensional example using the Cuntz algebra O_2 — a simple,
purely infinite C*-algebra with no traces and no finite-dimensional
representations. Key result: **BP forces type collapse** — constrained
extremal states must have type I GNS representations (finite-dimensional
commutant), despite generic states being type III factors.

New content in `discussion/barvinok_pataki_aou_proof.tex`:
- **Section 12** ("The Cuntz Algebra O_2"):
  - Definition, key properties (simple, purely infinite, no traces)
  - States as measures on Cantor space, Bernoulli (pure) vs Markov (type III)
  - Corollary: BP for Cuntz algebras (dim commutant ≤ m+1)
  - Type collapse table: m constraints → commutant structure
  - Truncation to M_{2^L} for numerical verification
  - Two numerical tables: gauge-invariant constraints, convergence with L
  - Markov constraints producing rank-2 mixed states at m=15

New Julia script `numerics/cuntz_sdp.jl`:
- 4 examples: gauge-invariant constraints, coherences, convergence, Markov vs Bernoulli
- Truncation levels L=2 to L=6 (d=4 to d=64)
- Rank stable across all truncation levels (confirms infinite-dim result)

### Numerical Results Tables (Sections 8–10)

Integrated SDP verification data into the conic inequality and spectrahedra sections:
- Table 5: Block-diagonal LMI spectrahedron (8 trials)
- Table 6: PPT entanglement witness ranks
- Table 7: Scalar inequality slack analysis (BP tight at m=15)
- Table 8: Compression constraint (Martin vs Lennart bounds)

### Email Follow-Up Integration (4 March 2026)

Integrated six threads from the post-meeting email discussion into the main
LaTeX document:

#### New Sections Added

- **Section 8: Cone-Valued Inequality Constraints** — Martin's extension:
  Theorem (BP for conic inequalities) proving d_F ≤ Σ dim(W_i) for constraints
  Φ_i(x) ≥_{W_i+} a_i. Remarks on equality vs inequality, simplex vs PSD cones,
  inactive constraint reduction (Timo), and slack variable formulation (Lennart).

- **Section 9: Application to Spectrahedra** — Block-diagonal LMI corollary
  (joint rank bound strictly tighter than individual), entanglement witness rank
  bound via PPT constraint (Martin's application suggestion).

- **Section 10: Inequality Constraints and the Role of Slack** — Lennart's
  PSD + scalar inequality example with slack variable analysis. Demonstrates that
  face reduction cannot control whether it reduces rank(X) or kills slacks.

- **Section 11: Obstructions to Nonlinear Extensions** — Lukas's manifold
  framework (differential version, open question on flow reaching lower face) and
  Alex's fundamental obstruction (interior shell counterexample showing linearity
  is essential).

#### Numerical Spectrahedra Examples

New Julia script `numerics/spectrahedra_bp.jl` with five examples:

| Example | Key Finding |
|---------|-------------|
| Block-diagonal LMI (n=10, two 4×4 blocks) | Joint bound r_F+r_G ≤ 4 tighter than individual r_F,r_G ≤ 4 |
| Lennart's scalar inequalities (PSD_6) | BP tight at m=15: face_dim=16=16. Timo's inactive reduction confirmed |
| PPT witness (2×2 system) | rank(W^Γ) bounded by m; rank(W) can be full while rank(W^Γ) stays low |
| Compression constraint (8×8, 3×3 subblock) | Martin and Lennart bounds both verified; slack approach tighter |

#### Email Discussion Summary

Written up in `discussion/email_extensions.md` with full analysis of all six
email threads (Martin, Timo, Lukas, Alex, Lennart) and connections to open questions.

---

## Work Completed (Session 3)

### Authorship and Credit

- **Timo Ziegler moved to first author** in `discussion/barvinok_pataki_aou_proof.tex`,
  with explicit credit in acknowledgements for the proof idea (perturbation–extremality
  argument) that was recognised to generalise to pointed cones and AOU spaces.
- Martin Plavala credited for recognising the general cone extension during the discussion.
- **Henrik Wilming added as author** (AG Osborne, ITP Hannover) and as participant
  in `discussion/meeting_summary.md`.

### Infinite-Dimensional Extension via Weis–Shirokov

Following an email from Alexander Stottmeister, incorporated **Theorem 3.3** from
Weis–Shirokov (J. Convex Anal. 28(3), 2021; arXiv:2003.14302) into the LaTeX write-up.

Key insight: the finite-dimensionality hypothesis on $d_F$ in the main theorem is
**automatically satisfied** for extreme points under finitely many constraints, by a
purely convex-geometric result that requires only a real vector space and a convex set
— no operator-algebraic notions needed. Operator algebras enter only to *interpret*
the finite-dimensional face as a commutant.

New content in `discussion/barvinok_pataki_aou_proof.tex`:
- **Section 5** ("The Infinite-Dimensional Case and the Role of $d_F$"):
  Theorem 8 (Weis–Shirokov), Corollary 9 (finiteness of $d_F$),
  Remark 10 (no operator-algebraic assumptions needed),
  Corollary 11 (infinite-dimensional Barvinok–Pataki without a priori finiteness of $d_F$)
- Updated abstract to reflect the infinite-dimensional extension
- Fixed Weis–Shirokov bibliography (was J. Math. Phys., corrected to J. Convex Anal.)
- Downloaded paper to `discussion/shirokov_weis_2021.pdf`

### Hydrogen Atom Example (Section 7)

Full worked physical example: hydrogen in a magnetic field with angular momentum
constraints (Jz, J², L²). Includes complete notation (Hilbert space, quantum numbers,
Coulomb + Zeeman Hamiltonians), constrained optimisation setup, rank bound tables for
both infinite-dimensional and finite-shell cases, and physical interpretation.

### Numerical SDP Verification

Julia program `numerics/hydrogen_sdp.jl` using JuMP + Mosek:
- Two-phase solver (optimal energy → degeneracy-breaking perturbation)
- Experiments across shells n=2,3,4,5 (dim 8–50), B ∈ [0.05, 5], m=0..5 constraints
- **BP bound verified**: all extreme-point optimizers satisfy r² ≤ m+1
- **Tight bound demonstrated**: rank-2 optimizers found when constraint values lie
  in the interior of the joint numerical range (no pure state feasible)
  - n=2: ⟨Jz⟩=0, ⟨J²⟩=2, ⟨L²⟩=1 → rank=2, eigenvalues (0.569, 0.431)
  - n=3: ⟨Jz⟩=0, ⟨J²⟩=4, ⟨L²⟩=3 → rank=2, eigenvalues (0.675, 0.325)
- Analytic-centre phenomenon documented (M=1/2 case)

5 TikZ/pgfplots figures integrated into the main paper:
1. Rank vs m bar chart (extreme point vs analytic centre)
2. Optimal energy vs B (n=3 shell)
3. Rank vs Hilbert space dimension (BP bound vs naïve r≤d)
4. Rank-2 eigenvalue spectrum (tight bound)
5. Table of genuine rank-2 optimizers

### Discussion Directory

- `discussion/barvinok_pataki_aou_proof.tex` — 10-page paper (compiles cleanly)
- `discussion/meeting_summary.md` — meeting notes from 3 March 2026
- `discussion/shirokov_weis_2021.pdf` — Weis–Shirokov paper
- `numerics/hydrogen_sdp.jl` — SDP solver + experiments
- `numerics/find_rank2.jl` — rank-2 optimizer search
- `numerics/plots.tex` — standalone plot compilation
- `numerics/data/` — 29 data files
- Audio/video recordings and board snapshots from the session

---

## Work Completed (Session 2)

### All 4 Genuine Math Errors Fixed

| Node | Error | Fix |
|------|-------|-----|
| **1.9.2** | Face dimension spectrum for M₂(C)⊕M₂(C) is {0,1,3,4,7}, NOT {0,3,7} | Fixed in session 1; remaining challenges resolved |
| **1.6.6** | "Trivial involution" wrong — real PSD cone uses transpose a*=aᵀ. M_n(R)^{sa} is not an algebra. Missing quaternionic case. | Amended with correct involution, proper notation A=M_n(R), and full Hurwitz classification (R/C/H) |
| **1.10.4** | Marginal constraint count: dim²−1, not dim² (trace auto-satisfied). Qubits: m=3k, not 4k. | Amended with explicit basis decomposition, overlapping subsystem caveat, correct count |
| **1.10.3** | KMV state polynomials are nonlinear in φ; theorem only handles affine constraints | Scope restricted to linear sub-case. Structural novelty (Artin-Wedderburn beyond rank) clarified |

### All 148 Challenges Resolved

| Category | Count | Resolution |
|----------|-------|------------|
| Genuine math errors | 4 | Statement amended via `af amend` |
| Critical logical gaps | 6 | Proof steps added/expanded in amendments |
| Wrong inference type ("assumption" for derived claims) | ~40 | Corrected in statement text (af amend only modifies text, not metadata) |
| Missing dependency declarations | ~40 | Declared in statement text |
| Imprecise statements / notation | ~30 | Clarified |
| Positive verification notes | ~10 | Acknowledged |
| Minor completeness / precision | ~18 | Addressed |

### All Proof Sections Addressed

- **Core proof (1.1–1.4):** Face-commutant positivity proof cleaned, domain convexity proved, interior-of-cone argument for dimension formula, Weis-Shirokov converted to external citation, hypothesis verification chain completed, Artin-Wedderburn with complete reducibility argument
- **Abelian recovery (1.5):** Circularity broken via L^∞ characterisation, GNS computation made explicit
- **Matrix recovery (1.6):** Involution corrected, quaternionic case added, extremality requirement explicit
- **Block-diagonal (1.7):** Inequivalence via central projections, cross-block novelty spelled out
- **Infinite-dim (1.8):** Normal state restriction justified rigorously, unbounded Hamiltonians via lsc generalised affine maps, SDP solver claim corrected (analytic centre, not extreme point)
- **Face gap (1.9):** Three-way case analysis for N=1,2,≥3; corrected spectrum table
- **Applications (1.10):** NPA finite-level vs convergence separated, KMV scope restricted, marginal count corrected, SDP solver misconception fixed
- **Type I extension (1.8.5):** Von Neumann algebra examples corrected to C*-algebras (C*_r(F_n))

### LaTeX Write-Up Produced

- `latex/main.tex` — 11-page self-contained document following af-tests style
- `latex/proof-tree-export.tex` — 78KB raw proof tree export
- Clearly marked as **not yet validated** with status box
- Includes: main theorem, 3-step proof, all recoveries, applications, error report, methodology, appendix with full proof tree

---

## Key Learnings

### 1. Proof Architecture

**The core proof is remarkably clean.** Three steps:
1. Face-commutant correspondence: dim(F) = dim_R(commutant^{sa}) - 1
2. Weis-Shirokov: dim(F) ≤ m for extreme points under m constraints
3. Combine: dim_R(commutant^{sa}) ≤ m + 1

All three steps are mathematically correct and were validated. The
difficulty was in supporting infrastructure.

### 2. Systematic Issues

**Every challenged node** had two metadata problems:
- Inference type "assumption" when it should have been a deduction
- Missing dependency declarations

These are limitations of the initial proof construction, not of the
mathematics. The `af amend` command only modifies statement text, so
the metadata fields remain as "assumption" — but the statement text
now explicitly identifies each node's logical role and dependencies.

### 3. Application Nodes Are Dangerous

The most serious errors (KMV scope mismatch, marginal constraint count,
SDP solver misconception, NPA conflation of finite/infinite levels)
were all in the applications section (1.10.x). These tend to have
hand-wavy connections to the main theorem.

### 4. Verification Process

- One verifier per node with strict isolation prevents cross-contamination
- Opus-level verifiers caught errors that would survive casual review
- Breadth-first verification surfaces structural issues early
- The adversarial protocol was effective: 4 genuine errors and 6 critical
  gaps found across 64 nodes

---

## Next Steps (Priority Order)

1. **Run verifier pass** on 49 pending nodes (`af jobs` → claim verifier
   jobs → `af accept`)
2. **Re-examine** any nodes that receive new challenges during re-verification
3. **Polish LaTeX** into a submission-ready paper
4. **Consider Lean 4 formalisation** of the core 3-step proof (leveraging
   the af-tests GNS infrastructure)

---

## Statistics

| Metric | Session 1 | Session 2 | Change |
|--------|-----------|-----------|--------|
| Nodes validated | 15 | 15 | — |
| Open challenges | 148 | **0** | −148 |
| Quality score | 39.4 | **69.4** | +30.0 |
| Genuine errors found | 4 | 4 | — |
| Amendments made | 3 | **38** | +35 |
| LaTeX pages | 0 | **11** | +11 |

- Definitions: 32
- External references: 14
- Proof nodes: 64 (11 Level-1, 53 Level-2)
- Verifier runs: 26 independent opus-level agents (session 1)
- Prover repair agents: 5 parallel + manual (session 2)
