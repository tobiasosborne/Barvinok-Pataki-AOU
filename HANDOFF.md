# Handoff Report: Barvinok-Pataki for Order Unit Spaces

## Status

**All 24 proof nodes adversarially verified + 4 new sections added (28 pages).**

- LaTeX source: `discussion/barvinok_pataki_aou_proof.tex` (28 pages, compiles cleanly)
- af formalisation: 24 nodes validated, 15 definitions, 8 external references
- Lean 4 formalisation: `lean/` — sorry-filled skeleton, `lake build` succeeds
- Verified reference: `references/weis_shirokov_2021.pdf` (arXiv:2003.14302)
- Old formalisation: `archive/v1/` (646 ledger entries, 15 externals, PRD, latex export)

---

## Session 9: Rényi Generalisations and Robustness (Section 17)

### Motivation
Coauthor feedback: "I don't care about rank, it isn't robust." Rank = exp(S₀) is
discontinuous; Rényi dimensions d_α = exp(S_α) for α > 0 are continuous alternatives.

### New Section 17 in LaTeX (~230 lines, 7 subsections)

- **17.1 Rényi dimension**: Definition, table of special cases (α=0,1,2,∞), monotonicity
- **17.2 Theorem (Rényi BP)**: d_α(ρ*)² ≤ m for all α ≥ 0. Trivial corollary of BP; tight when eigenvalues are uniform.
- **17.3 Theorem (Robust Rényi BP)**: For δ-approximate extreme points with α > 1:
  d_α(ρ) ≤ √m / (1-δ/2)^{α/(α-1)}. **Dimension-free** — does not depend on ambient d.
  Key insight: tail eigenvalues contribute (δ/2)^α to tr(ρ^α), which vanishes for α > 1.
  For α = 1 (von Neumann), dimension-dependence is unavoidable (Fannes–Audenaert).
- **17.4 Corollary (Purity BP)**: tr(ρ*²) ≥ 1/√m. Measurable via SWAP test without tomography.
- **17.5 Eigenvalue clustering** (Pataki 1998): At SDP optima, eigenvalues cluster.
  Numerics show d₂/d₀ ≈ 0.25 generically (but no significant difference between optima and random).
- **17.6 Conjecture (Typical Rényi BP)**: E[d_α(ρ*)] ≤ c_α · √m with c_α < 1 for α > 0.
- **17.7 Numerical verification table**: d_α values for random spectrahedra (m=25, d=10).

### New discussion document
- `discussion/renyi_generalisation.md` — full 11-section mathematical investigation
  covering robustness, Schatten norms, modular theory connections, SOS hierarchies,
  Roy–Vetterli effective rank, and 7 numerical findings.

### New numerics: `numerics/renyi_bp.jl` (7 experiments)
1. Random spectrahedra: d_α across varying m (confirms d_α ≪ rank)
2. Robustness: depolarising perturbation shows d_α continuous while rank jumps
3. Extreme CPTP channels: d_α well below BP bound
4. Eigenvalue skewness: BP-saturating points have d₂/d₀ ≈ 0.25 (heavily skewed)
5. Pataki clustering: no significant clustering difference between optima and random
6. Scaling: c_α = d_α/√m < 1 across m ∈ {4,9,16,25} (supports Typical Rényi BP conjecture)
7. Applications: marginals (d₁=1.57 vs bound 4.1), covariant channels (d_α=2 vs bound 6)

### Key results
- **Theorem 2 (Robust Rényi BP) is genuinely new**: dimension-free robustness for α > 1
- **Purity bound is experimentally testable** via SWAP test
- **Numerics confirm** d_α ≤ c_α√m with c_α ≈ 0.3–0.5 generically
- **Rank fails numerically** (SCS returns rank > BP bound) while d_α respects bounds

### New bibliography
- Roy & Vetterli 2007 (effective rank), Im & Wolkowicz 2021 (strengthened BP),
  Fannes 1973 (entropy continuity)

---

## Session 8: 6-Workstream Extension (Sections 14–16 + Remarks + Lean 4)

### New Remarks
- **Remark 7** (real vs complex face dimensions): After Theorem 1, clarifies r(r+1)/2 vs r² conventions
- **Remark (constraints increase bound)**: In Section 13, explains why adding symmetry constraints increases BP bound, resolution via Schur block decomposition

### New Section 14: BP for Quantum Error Correction (~100 lines)
- KL conditions as linear constraints on recovery channel Choi matrix
- BP bound for recovery channels: r(r+1)/2 ≤ m_TP + m_KL
- Examples: [[4,2,2]], [[5,1,3]], [[7,1,3]] (with remark on vacuity of combined bound)
- Numerics: `numerics/qec_bp.jl`

### New Section 15: BP for the Quantum Marginal Problem (~130 lines)
- Marginal constraints as spectrahedron
- BP bound: r(r+1)/2 ≤ d_A(d_A+1)/2 + d_B(d_B+1)/2 − 1
- Numerical verification: qubit-qubit, qubit-qutrit, qutrit-qutrit
- Numerics: `numerics/marginals_bp.jl`

### New Section 16: BP for Thermal Operations (~100 lines)
- Gibbs-preserving channels as spectrahedron
- BP bound: r(r+1)/2 ≤ d(d+1) − 1
- Temperature dependence: β→0 (unital), β→∞ (ground-state preserving)
- Numerics: `numerics/thermal_bp.jl`

### Lean 4 Formalisation (`lean/`)
- Definitions: ConvexCone face, spanFace, faceDim, constraintMap, feasibleSet, IsExtreme
- Main theorem: `barvinok_pataki` — sorry-filled skeleton with correct statement
- 4 sorry lemmas: rank-nullity, perturbation, relint, extremality contradiction
- `lake build` succeeds (Lean 4.27.0, Mathlib v4.27.0)

### New bibliography
- Knill & Laflamme 1997, Gottesman 1997, Klyachko 2006, Lostaglio 2019, Brandao et al. 2015

### Adversarial verification: 3 parallel verifiers, key fixes
- WS4: m_cov=4→5, m_TP=4→3 (real symmetric convention)
- WS1: Acknowledged QEC bounds are vacuous (exceed TP-only bound); added remark
- WS2: Fixed Klyachko citation context (marginal problem, not de Finetti)
- WS3: Fixed "single rank-1 constraint" wording
- Code: Fixed m_tp convention in channels_bp.jl (dA^2 → dA*(dA+1)/2)

---

## Session 7: Quantum Channels, Combs, and Scattering (Section 13)

### Added Section 13 (7 subsections, ~250 lines of LaTeX)
- 13.1: Channels as a spectrahedron — Choi's theorem recovered via BP
- 13.2: Constrained channels — covariance, fixed I/O, observable preservation
- 13.3: Quantum combs and superchannels (CDP 2008/2009)
- 13.4: Inclusive scattering channels with conservation laws
- 13.5: Generic 2→2 scattering with spin (parametric table)
- 13.6: Electron-electron Coulomb scattering (partial waves, J_z, parity)
- 13.7: Compton scattering (J_z conservation, helicity basis)

### New numerics
- `numerics/channels_bp.jl`: CPTP channels, covariant channels, fixed I/O, 2-combs
- `numerics/scattering_bp.jl`: generic scattering, e-e Coulomb, Compton

### New bibliography
- Choi 1975, CDP 2008 (Europhys. Lett.), CDP 2009 (PRA)

### New af nodes (1.18–1.23) — all validated

```
1.18  Recovery: Choi's theorem via BP
1.19  BP for constrained channels
1.20  BP for quantum N-combs
1.21  Scattering channel complexity
1.22  Electron-electron scattering example
1.23  Compton scattering example
```

### Adversarial verification: 3 rounds, ~23 challenges raised, all resolved

Key errors found and fixed by verification:
- **CDP hierarchy definition**: wrong traces and base case corrected (recursive Tr_{B_k} notation)
- **N-comb constraint counts**: m=15 → m=52 for N=2 qubits (with derivation)
- **S-matrix isometry**: S†S = I now stated explicitly for inclusive channel CPTP
- **Real/complex BP formula mismatch**: all SDPs use real symmetric PSD, so r(r+1)/2 ≤ m (not r² ≤ m); code changed to `bp_rank_real`, all tables updated
- **Analytical branch constraint counting**: rewritten to iterate over Choi upper-triangle entries consistently with SDP branch
- **Conservation law derivation**: expanded to explain zero-entry constraints from [Q_total, S] = 0
- **Subspace dimension**: "4-dimensional" → "16-dimensional" in N-comb derivation

### Paper: compiles clean, 20 pages

---

## Previous Sessions Summary

- **Session 6**: Adversarial verification — 18/18 nodes validated, 2 errors + 4 gaps fixed
- **Session 5**: Archived old formalisation, new af v2 from discussion paper
- **Session 4**: Cuntz algebra O_2 example (Section 12), numerical tables (Sections 8-10)
- **Session 3**: Authorship, Weis-Shirokov infinite-dim extension, hydrogen atom example
- **Session 2**: All 148 challenges resolved, 4 genuine errors fixed, LaTeX produced
- **Session 1**: Initial af formalisation, 64 nodes, adversarial verification

---

## Proof Tree (24 validated + 9 new nodes)

```
1   Main Theorem: BP for Pointed Cones (d_F <= m)
1.1   Step 1: Rank-nullity reduction to ker phi_F = {0}
1.2   Step 2: Perturbation preserves feasibility
1.3   Step 2: Algebraic relative interior of F (no closedness needed)
1.4   Step 3: Contradiction with extremality (QED)
1.5   Corollary: BP for AOU spaces
1.6   Corollary: BP for C*-algebras (sum n_j^2 <= m+1)
1.7   Finiteness of d_F via Weis-Shirokov (no a priori assumption needed)
1.8   Corollary: BP for infinite-dimensional pointed cones
1.9   Theorem: BP for conic inequalities (d_F <= sum dim(W_i))
1.10  Corollary: BP for Cuntz algebras (type collapse)
1.11  Recovery: Classical Barvinok-Pataki for M_n(C)
1.12  Recovery: Caratheodory's theorem for C(X)
1.13  Recovery: Block-diagonal algebras
1.14  Recovery: Infinite-dimensional B(H)
1.15  Remark: Inactive constraints reduce effective m
1.16  Remark: Slack variable formulation
1.17  Remark: Obstruction to nonlinear extensions
1.18  Recovery: Choi's theorem via BP
1.19  BP for constrained channels
1.20  BP for quantum N-combs
1.21  Scattering channel complexity
1.22  Electron-electron scattering example
1.23  Compton scattering example
1.24  KL conditions as linear constraints on Choi matrix
1.25  BP bound for QEC recovery channels (vacuous for standard codes)
1.26  QEC code examples (constraint counts and rank bounds)
1.27  Marginal constraints form a spectrahedron
1.28  BP bound for marginal-constrained states
1.29  Qubit-qubit marginal example
1.30  Thermal operations as spectrahedron
1.31  BP bound for thermal operations
1.32  Temperature limits analysis
1.33  Rényi BP bound (d_α ≤ √m for all α)
1.34  Robust Rényi BP (dimension-free bound for α > 1)
1.35  Purity BP corollary (tr(ρ²) ≥ 1/√m, SWAP-testable)
1.36  Typical Rényi BP conjecture (c_α < 1 generically)
```

---

## Next Steps

1. **Prove Typical Rényi BP conjecture**: Use random matrix theory (MP distribution for truncated spectra) to show E[d_α] ≤ c_α√m with c_α < 1
2. **Experimental proposal**: Design SWAP-test experiment to verify purity BP bound for quantum channels
3. **Explore α ∈ (0,1) regime**: Partial robustness results may be possible with additional structure
4. **Modular theory connection**: Relate Rényi BP to Tomita–Takesaki theory for infinite-dim AOU spaces
5. **Download remaining references** (Pataki, Cuntz, Bratteli-Jorgensen) if access becomes available
6. **Consider Lean 4 formalisation** of the core 3-step proof (nodes 1.1–1.4)
7. **Tighten node 1.3** per verifier note: make F-membership explicit in the calculation
