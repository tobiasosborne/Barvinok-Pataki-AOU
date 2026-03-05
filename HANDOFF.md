# Handoff Report: Barvinok-Pataki for Order Unit Spaces

## Status

**All 18 proof nodes adversarially verified (100% validated, 0 open challenges).**

- LaTeX source: `discussion/barvinok_pataki_aou_proof.tex` (17 pages, compiles cleanly)
- af formalisation: 18 nodes validated, 15 definitions, 8 external references (2 verified via string match)
- Ledger: 129 entries (prover amendments + 2 rounds of adversarial verification)
- Verified reference: `references/weis_shirokov_2021.pdf` (arXiv:2003.14302)
- Old formalisation: `archive/v1/` (646 ledger entries, 15 externals, PRD, latex export)

---

## Work Completed (Session 5)

### Archived Old Formalisation

Moved the original af formalisation (sessions 1-2) to `archive/v1/`:
- `bp-order-unit-prd.md` (PRD)
- `ledger/` (646 entries)
- `externals/` (15 references)
- `meta.json`
- `latex/` (main.tex, proof-tree-export.tex)

### New af Formalisation (v2)

Created from scratch based on `discussion/barvinok_pataki_aou_proof.tex`.

#### Definitions (15)
ordered_vector_space, pointed_cone, archimedean_order_unit_space,
face_of_cone, state_space, GNS_construction, face_commutant_correspondence,
extreme_point, spectrahedron, cuntz_algebra, generalised_affine_map, rank_nullity

#### External References (8)
| Reference | Status |
|-----------|--------|
| Weis-Shirokov Theorem 2 (arXiv:2003.14302) | VERIFIED (local PDF, string matched) |
| Weis-Shirokov Corollary 5 (arXiv:2003.14302) | VERIFIED (local PDF, string matched) |
| Alfsen-Shultz Face-Commutant Correspondence (2003) | Book, not downloadable |
| Pataki SDP Rank Bound (1998) | Paywalled journal |
| Barvinok Convexity (2002) | Book, not downloadable |
| Rockafellar Convex Analysis (1970) | Book, not downloadable |
| Cuntz Simple C*-algebras (1977) | Paywalled (Springer) |
| Bratteli-Jorgensen IFS and Cuntz (1999) | Download failed (old arXiv format) |

#### Proof Tree (18 nodes)

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
```

### Key Design Decisions

1. **Conjecture = Theorem 2 (general pointed cones)**, not the C*-algebra version.
   The C*-algebra corollary (node 1.6) depends on the general theorem plus the
   Alfsen-Shultz external reference.

2. **Weis-Shirokov is the only proof-critical external reference** (used in node 1.7
   for infinite-dimensional finiteness of d_F). All other externals are background/motivation.

3. **Strict verification rule**: references downloaded to `references/` and results
   string-matched via `pdftotext`. Only Weis-Shirokov passes this bar; others are
   books or paywalled journals marked accordingly.

---

## Session 6: Adversarial Verification (4 parallel verifiers)

### Findings: 2 errors, 4 gaps, 5 warnings — all fixed

**Errors fixed:**
- E1: Off-by-one in Thm 2 "equivalently" (dim(F)≤m−1 → dim(F)≤m) + Corollary 4
- E2: "Archimedean implies pointed" false; added "pointed" to AOU definition

**Gaps fixed:**
- G1: Closedness gap in Steps 2–3. Entire proof replaced with direct relint perturbation (3 steps, no closedness needed)
- G2: "x* on proper face" unjustified — eliminated by simplified proof
- G3: d_F finiteness missing for inf-dim C*-algebras — added forward ref to Weis-Shirokov
- G4: Inactive constraint "small |t|" argument insufficient — now cites slack formulation

**Warnings addressed:**
- Cuntz table m=3 row expanded (was incomplete)
- B(H) recovery now cites trace-class route directly
- "Polyhedron" → "cone" in nonlinear obstruction
- Section 3 preamble updated to match simplified proof

**Paper:** compiles clean, 17 pages

### Round 2: Re-verification (4 parallel verifiers via af CLI)

All 18 nodes re-verified after fixes. Results:
- **18/18 accepted** (0 challenges remaining)
- 1 cosmetic challenge raised (node 1.9 stale step numbering) — resolved and re-accepted
- Node 1.3 accepted with note: explicit calculation shows ∈ V⁺, membership in F relies on abstract relint principle (correct but could be more transparent)
- Node 1.8 accepted with note: Thm 2 redundant since Weis-Shirokov alone gives d_F ≤ m (not an error)

**Ledger:** 129 entries (000043–000129: amendments, claims, verifications, acceptances)

---

## Previous Sessions Summary

- **Session 4**: Cuntz algebra O_2 example (Section 12), numerical tables (Sections 8-10)
- **Session 3**: Authorship, Weis-Shirokov infinite-dim extension, hydrogen atom example
- **Session 2**: All 148 challenges resolved, 4 genuine errors fixed, LaTeX produced
- **Session 1**: Initial af formalisation, 64 nodes, adversarial verification

---

## Next Steps

1. **Download remaining references** (Pataki, Cuntz, Bratteli-Jorgensen) if access becomes available
2. **Add sub-proofs** for corollary nodes (1.5, 1.6, 1.9, 1.10) if deeper formalisation desired
3. **Consider Lean 4 formalisation** of the core 3-step proof (nodes 1.1–1.4)
4. **Tighten node 1.3** per verifier note: make F-membership (not just V⁺-membership) explicit in the calculation
