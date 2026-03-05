/-
  BPForCones/MainTheorem.lean — The Barvinok–Pataki theorem for pointed cones

  Main result: If x is an extreme point of the feasible set
    C = { x ∈ V_pos | α_j(x) = b_j, j = 1..m }
  then d_F = dim(span(face(x))) ≤ m.

  A "pointed cone" in the paper corresponds to Mathlib's
  `ConvexCone.Salient`: ∀ x ∈ C, x ≠ 0 → -x ∉ C.
-/
import BPForCones.Basic
import Mathlib.LinearAlgebra.Dimension.Finrank
import Mathlib.LinearAlgebra.FiniteDimensional.Defs

set_option autoImplicit true

open Submodule BPCone

namespace BPCone

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-! ## Step 1: Rank-nullity for the constraint map -/

/-- If the constraint map φ_F has trivial kernel, then d_F ≤ m. -/
theorem faceDim_le_of_ker_trivial
    (C : ConvexCone ℝ V) (x : V)
    (α : Fin m → (V →ₗ[ℝ] ℝ)) (b : Fin m → ℝ)
    (hfin : FiniteDimensional ℝ (spanFace C x))
    (hker : LinearMap.ker (constraintMap C x α) = ⊥) :
    faceDim C x ≤ m := by
  sorry

/-! ## Step 2: Perturbation preserves feasibility -/

/-- If D ∈ ker φ_F, then x + tD satisfies all constraints. -/
theorem perturbation_preserves_constraints
    (C : ConvexCone ℝ V) (x : V)
    (α : Fin m → (V →ₗ[ℝ] ℝ)) (b : Fin m → ℝ)
    (hx : x ∈ feasibleSet C α b)
    (D : spanFace C x)
    (hD : D ∈ LinearMap.ker (constraintMap C x α))
    (t : ℝ) :
    ∀ (j : Fin m), α j (x + t • (D : V)) = b j := by
  sorry

/-- If x ∈ relint(face(x)) and D ∈ span(face(x)), then ∃ ε > 0 with x ± εD ∈ face. -/
theorem relint_perturbation_in_face
    (C : ConvexCone ℝ V) (x : V)
    (hx : x ∈ C)
    (D : spanFace C x)
    (hD : D ≠ 0) :
    ∃ ε : ℝ, 0 < ε ∧
      (x + ε • (D : V)) ∈ face C x ∧
      (x - ε • (D : V)) ∈ face C x := by
  sorry

/-! ## Step 3: Contradiction with extremality -/

/-- Extremality implies ker φ_F = {0}. -/
theorem ker_trivial_of_extreme
    (C : ConvexCone ℝ V) (x : V)
    (α : Fin m → (V →ₗ[ℝ] ℝ)) (b : Fin m → ℝ)
    (hsal : C.Salient)
    (hfin : FiniteDimensional ℝ (spanFace C x))
    (hx : x ∈ feasibleSet C α b)
    (hext : IsExtreme (feasibleSet C α b) x) :
    LinearMap.ker (constraintMap C x α) = ⊥ := by
  sorry

/-! ## Main theorem -/

/-- **Barvinok–Pataki for pointed cones.**
    Let V be a real vector space, C ⊆ V a salient (pointed) convex cone, and
    α₁, ..., αₘ linear functionals. If x is an extreme point of
    { x ∈ C | α_j(x) = b_j }, then d_F ≤ m. -/
theorem barvinok_pataki
    (C : ConvexCone ℝ V) (x : V)
    (α : Fin m → (V →ₗ[ℝ] ℝ)) (b : Fin m → ℝ)
    (hsal : C.Salient)
    (hfin : FiniteDimensional ℝ (spanFace C x))
    (hx : x ∈ feasibleSet C α b)
    (hext : IsExtreme (feasibleSet C α b) x) :
    faceDim C x ≤ m := by
  exact faceDim_le_of_ker_trivial C x α b hfin
    (ker_trivial_of_extreme C x α b hsal hfin hx hext)

end BPCone
