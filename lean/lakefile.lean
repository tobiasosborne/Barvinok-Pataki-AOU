import Lake
open Lake DSL

package «bp-for-cones» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.27.0"

@[default_target]
lean_lib BPForCones where
  srcDir := "."
  roots := #[`BPForCones.Basic, `BPForCones.MainTheorem]
