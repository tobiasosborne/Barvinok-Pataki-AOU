# Barvinok-Pataki for Order Unit Spaces

Adversarial formalisation of an abstract Barvinok-Pataki theorem for state
spaces of Archimedean order unit spaces, unifying three previously
disconnected literatures:

1. **SDP optimisation** — Barvinok-Pataki rank bounds
2. **Abstract convex geometry** — Weis-Shirokov face-dimension theorem
3. **Operator algebras** — Alfsen-Shultz face-commutant correspondence

## Main Result

Let **A** be a C\*-algebra (or Archimedean order unit space) and let
a₁, …, aₘ ∈ A^sa be self-adjoint elements. If φ is an extreme point of
the constrained state space C = {φ ∈ S(A) : φ(aⱼ) = βⱼ}, then:

**Σ nⱼ² ≤ m + 1**

where π_φ ≅ ⊕ σⱼ^{⊕nⱼ} is the GNS decomposition into irreducible
representations with multiplicities nⱼ.

This simultaneously recovers:
- Carathéodory's theorem (A = C(X) abelian)
- Classical Barvinok-Pataki for complex SDP (A = Mₙ(ℂ))
- A new mixed bound for block-diagonal algebras (A = ⊕ M_{dᵢ}(ℂ))
- Weis-Shirokov's theorem on energy-constrained quantum states (A = B(H))

## Structure

This repository contains an **adversarial proof formalisation** built with
the [`af` CLI tool](https://github.com/af-framework/af), together with a
LaTeX write-up of the proof.

```
.
├── bp-order-unit-prd.md          # Full mathematical PRD (problem specification)
├── HANDOFF.md                    # Verification report with learnings
├── README.md                     # This file
├── latex/
│   ├── main.tex                  # LaTeX write-up (11 pages)
│   ├── main.pdf                  # Compiled PDF
│   └── proof-tree-export.tex     # Full af export of proof tree
├── ledger/                       # Append-only proof ledger (af framework)
├── meta.json                     # Proof metadata (af framework)
├── defs/                         # Definition files
├── externals/                    # External reference files
└── nodes/                        # Proof node data (af framework)
```

### Proof Tree (64 nodes)

```
1   Main Theorem (Barvinok-Pataki for Order Unit Spaces)
├── 1.1   Face-Commutant Correspondence (Alfsen-Shultz)
│   └── 1.1.1–1.1.6  Well-definedness, affinity, injectivity,
│                     surjectivity, dimension consequence
├── 1.2   Abstract Face-Dimension Bound (Weis-Shirokov)
│   └── 1.2.1–1.2.4  Single/multi-constraint, extremal specialisation
├── 1.3   Main Dimension Bound (3-step core proof)
│   └── 1.3.1–1.3.5  Hypothesis check, bridge, bound, conclusion, mixed
├── 1.4   Commutant Decomposition → Multiplicity Bound
│   └── 1.4.1–1.4.5  Finite-dim, Artin-Wedderburn, decomposition, dim calc
├── 1.5   Recovery: Abelian → Carathéodory
├── 1.6   Recovery: Matrix Algebra → Classical BP
├── 1.7   Recovery: Block-Diagonal → Mixed Bound
├── 1.8   Recovery: Infinite-Dim → Weis-Shirokov
├── 1.9   Face Dimension Gap & Automatic Purity
├── 1.10  Applications (NPA, state polynomials, quantum marginals)
└── 1.11  Consistency Verification (counterexample check)
```

### Definitions (32)

All key mathematical objects are defined at a level readable by a graduate
student in functional analysis, with precise references to Alfsen-Shultz,
Takesaki, Weis-Shirokov, Barvinok, and Pataki.

## Verification Status

| Metric | Value |
|--------|-------|
| Nodes verified | 53/53 leaf nodes (100%) |
| Nodes accepted | 15 (23%) |
| Open challenges | **0** (all 148 resolved) |
| Genuine math errors found | 4 (all corrected) |
| Critical gaps found | 6 (all filled) |
| Quality score | 69.4/100 |

All 148 challenges raised during adversarial verification have been
resolved. The 49 pending nodes await a second verifier pass for
formal acceptance.

### Errors Found and Corrected

1. **Face dimension spectrum** — M₂(ℂ)⊕M₂(ℂ) has spectrum
   {0,1,3,4,7}, not {0,3,7}. The gap {1,2} only holds for simple algebras.
2. **"Trivial involution"** — The real PSD cone uses the transpose
   involution (a\*=aᵀ), not the identity. Added the full Hurwitz
   classification (real/complex/quaternionic).
3. **KMV scope mismatch** — State polynomials are nonlinear in φ; the
   main theorem only handles affine constraints. Application restricted
   to linear sub-case.
4. **Marginal constraint count** — Each marginal gives dim²−1 constraints
   (not dim²); trace is auto-satisfied. For qubits: m = 3k, not 4k.

See [HANDOFF.md](HANDOFF.md) for the full verification report.

## LaTeX Write-Up

An 11-page write-up is available in [`latex/main.tex`](latex/main.tex)
(compiled PDF: [`latex/main.pdf`](latex/main.pdf)). It includes:

- Main theorem statement and 3-step proof
- Recovery of all four classical results
- Face dimension gap and automatic purity theorem
- Applications (NPA, quantum marginals)
- Documentation of errors found during verification
- Appendix with full proof tree

**Note:** The write-up is clearly marked as a formalisation report, not a
finished paper. The proof structure is complete but formal validation is
pending.

## How to Explore

```bash
# Install af CLI (Go)
go install github.com/af-framework/af@latest

# View proof status
af status

# Check proof health and metrics
af health
af metrics
af progress

# List definitions
af defs

# View a specific node
af get 1.3    # Main dimension bound

# See available work (49 verifier jobs)
af jobs

# View challenge history on a node
af get 1.9.2  # The corrected face spectrum

# Export proof tree
af export --format latex -o proof.tex
af export --format markdown
```

## References

- Barvinok (2001). *A remark on the rank of PSD matrices subject to affine constraints*
- Pataki (1998, 2000). *On the rank of extreme matrices in SDP*
- Weis & Shirokov (2021). *The face generated by a point, generalized affine constraints, and quantum theory*
- Alfsen & Shultz (2001, 2003). *State Spaces of Operator Algebras*
- Klep, Magron & Volčič (2023). *State polynomials*
- Navascués, Pironio & Acín (2008). *NPA hierarchy*

## License

Apache License 2.0
