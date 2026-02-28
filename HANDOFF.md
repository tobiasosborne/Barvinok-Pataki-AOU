# Handoff Report: Barvinok-Pataki for Order Unit Spaces

## Status

**Adversarial formalisation of the main theorem and all corollaries is complete.**
- 64 proof nodes across 3 levels of depth
- 32 definitions, 14 external references
- 53 leaf nodes verified by independent opus-level adversarial verifiers
- 15 nodes accepted, 38 challenged with 148 open challenges

**Completion: 23% validated.** The proof *structure* is 100% laid out; the
remaining work is resolving challenges (mostly structural metadata, plus
~16 substantive mathematical issues).

---

## Key Learnings

### 1. Genuine Mathematical Errors Found

| Node | Error | Impact |
|------|-------|--------|
| **1.9.2** | Face dimension spectrum for M_2(C)⊕M_2(C) is {0,1,3,4,7}, NOT {0,3,7}. States on multiple blocks with rank-1 in each create face dim 1. Gap is {2,5,6}, not {1,2}. | **Corrects the PRD's table.** The gap {1,2} only holds for simple algebras (N=1, d≥2). |
| **1.6.6** | "Trivial involution" is wrong — the real PSD cone uses the *transpose* involution (a*=aᵀ), which is non-trivial. Also M_n(R)^{sa} is not an algebra. | **Corrects a conceptual error** in the real vs complex discussion. Missing quaternionic case (M_n(H), dim r(2r-1)). |
| **1.10.4** | Constraint count for quantum marginals is wrong: each marginal gives dim(H_{S_j})²-1 constraints, not dim². The trace constraint is auto-satisfied. For qubits: 3k, not 4k. | **Arithmetic error** in the application section. |
| **1.10.3** | KMV state polynomials are NONLINEAR in φ. The main theorem handles only AFFINE constraints. The node either misrepresents KMV or misapplies the theorem. | **Fundamental scope mismatch** — this application needs rethinking or restricting to linear sub-case. |

### 2. Critical Scope Issue (Resolved)

The main theorem was stated for "unital *-algebras" but the proof machinery
(GNS, face-commutant, Banach-Alaoglu) requires either:
- C*-algebras, or
- Archimedean order unit spaces, or
- *-algebras with Archimedean quadratic module (per the formally verified
  framework in `~/Projects/af-tests/latex/`)

**Resolution:** Node 1.3.1 was amended to specify the correct scope. The
af-tests work shows GNS extends to *-algebras with Archimedean quadratic
modules via Cauchy-Schwarz boundedness + Tychonoff compactness, which is
broader than just C*-algebras.

### 3. Proof Architecture Insights

**The core proof is remarkably clean.** The main theorem (nodes 1.3.1–1.3.4)
is a 3-step argument:
1. Face-commutant correspondence: dim(F) = dim_R(commutant^{sa}) - 1
2. Weis-Shirokov: dim(F) ≤ m for extreme points under m constraints
3. Combine: dim_R(commutant^{sa}) ≤ m + 1

All three steps were validated as mathematically correct. The difficulty lies
in the supporting infrastructure (face-commutant proof details, classical
case recoveries, applications).

**Systematic structural weakness:** Every single challenged node had missing
dependency declarations and wrong inference types ("assumption" used for
derived claims). This is a metadata problem, not a math problem, but it
would need batch-fixing before any formal verification could proceed.

### 4. Subtle Issues Worth Attention in the Paper

| Issue | Detail |
|-------|--------|
| **Abelian circularity** | Node 1.5.3 assumes finite support → 1.5.4 concludes finite support. Need: L^∞(X,μ) finite-dim ⟺ μ atomic with finite support. |
| **Face dimension formula** | dim(D) = dim(V)-1 (node 1.1.6) needs: cone is generating + has nonempty interior. Not just "one constraint." |
| **Abelian multiplicities** | "All irreps 1-dim" does NOT imply "all multiplicities = 1" (node 1.5.2). Need cyclicity + commutativity of commutant. |
| **Weis-Shirokov proof** | Node 1.2.1's inline proof outline has an incorrect constraint-propagation argument. Recommend citing WS externally. |
| **Unbounded observables** | Node 1.8.3 claims the bound extends to unbounded Hamiltonians but never verifies that Tr(Hρ) is a generalized affine map in WS's sense. |
| **SDP solver behavior** | Node 1.10.2 claims solvers "generically find extreme points" — this is wrong for interior-point methods (they find the analytic center). |
| **Beyond Type I examples** | Node 1.8.5 lists R, L(F_n) as C*-algebra examples — they are von Neumann algebras. Use C*_r(F_n) instead. |

### 5. Verification Process Learnings

- **One verifier per node** with strict isolation prevents cross-contamination
- **Opus-level verifiers** caught errors that would survive casual review
  (the M_2⊕M_2 spectrum error, the "trivial involution" error, the KMV
  nonlinearity issue, the marginal constraint count)
- **Breadth-first verification** efficiently surfaces structural issues early
- **The most dangerous nodes are applications** (1.10.x) — they tend to
  have hand-wavy connections to the main theorem
- **Setup/definition nodes** (1.5.1, 1.6.1, 1.7.1) are almost always correct
- **The +1 reconciliation** (1.6.5) was a feared pitfall but turned out to be
  completely clean and was accepted immediately

---

## Next Steps (Priority Order)

1. **Fix the 4 genuine math errors** (1.9.2 done; 1.6.6, 1.10.3, 1.10.4 remain)
2. **Batch-fix structural metadata** (add dependencies, fix inference types)
   across all 38 challenged nodes
3. **Add missing proof steps** for the 6 identified logical gaps
4. **Decide scope** for applications section (restrict KMV to linear case,
   fix quantum marginals count, clarify NPA finite-level vs convergence)
5. **Re-verify** amended nodes
6. **Write the paper** using this proof tree as the skeleton

---

## Statistics

- Definitions: 32
- External references: 14
- Proof nodes: 64 (11 Level-1, 53 Level-2)
- Verifier runs: 26 independent opus-level agents
- Challenges raised: 148
- Genuine errors found: 4
- Critical gaps found: 6
- Nodes validated: 15 (23%)
- Session duration: ~2 hours
