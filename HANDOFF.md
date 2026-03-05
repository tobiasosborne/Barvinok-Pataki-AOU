# Handoff Report: Barvinok-Pataki for Order Unit Spaces

## Status

**New af formalisation (v2) based on the 17-page paper. Old formalisation archived to `archive/v1/`.**

- LaTeX source: `discussion/barvinok_pataki_aou_proof.tex` (17 pages, compiles cleanly)
- af formalisation: 18 nodes, 15 definitions, 8 external references (2 verified via string match)
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
1.3   Step 3: Pointedness forces exit from cone
1.4   Step 4: Contradiction with extremality (QED)
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

## Previous Sessions Summary

- **Session 4**: Cuntz algebra O_2 example (Section 12), numerical tables (Sections 8-10)
- **Session 3**: Authorship, Weis-Shirokov infinite-dim extension, hydrogen atom example
- **Session 2**: All 148 challenges resolved, 4 genuine errors fixed, LaTeX produced
- **Session 1**: Initial af formalisation, 64 nodes, adversarial verification

---

## Next Steps

1. **Run verifier pass** on the 18 nodes (`af jobs` -> claim verifier jobs -> `af accept`)
2. **Download remaining references** (Pataki, Cuntz, Bratteli-Jorgensen) if access becomes available
3. **Add sub-proofs** for corollary nodes (1.5, 1.6, 1.9, 1.10) if deeper formalisation desired
4. **Consider Lean 4 formalisation** of the core 4-step proof (nodes 1.1-1.4)
