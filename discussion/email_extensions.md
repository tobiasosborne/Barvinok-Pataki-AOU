# Email Follow-Up: Extensions of the Barvinok-Pataki Theorem

**Date:** 4 March 2026 (day after blackboard discussion)
**Participants:** Martin, Timo, Lukas, Alex, Lennart, Tobias

## Summary of Proposed Extensions

### 1. Martin: Cone-Valued Inequality Constraints

**Core idea:** Extend from scalar-valued constraints f(x) >= a to conic inequality
constraints Phi_i: V -> W_i with Phi_i(x) >=_{W_i+} a_i.

**Setup:** Let (V, V+) be an ordered vector space with pointed cone. For a finite
index set I, let (W_i, W_i+) be ordered vector spaces with linear maps
Phi_i: V -> W_i. Impose constraints Phi_i(x) >= a_i (w.r.t. W_i+).

**Claimed bound:** d_F <= m where m = sum_i dim(W_i).

**Proof sketch:** Same perturbation-extremality argument. Find D in ker(Phi|_{span(F)})
where Phi: span(F) -> bigoplus_i W_i. Since Phi(D) = 0, the constraint
Phi(x + tD) = Phi(x) >= a is preserved. Walk to the boundary of V+ in the
direction of D. By rank-nullity and extremality: d_F <= sum dim(W_i).

**Key distinction for inequalities:**
- **Equalities:** Phi_i(x) = a_i gives dim(W_i) scalar constraints. Equivalent to
  the scalar formulation.
- **Inequalities with simplex cones** (R_+^k): Equivalent to k scalar inequality
  constraints. The conic bound recovers the scalar bound exactly.
- **Inequalities with non-simplex cones** (e.g., PSD_k): The conic bound d_F <= dim(W_i)
  is potentially tighter than converting to scalar constraints, because the PSD structure
  couples the constraints.

**Application:** Upper bounding the rank of entanglement witnesses. If the witness
must satisfy a PPT constraint (W^Gamma >= 0), this is a conic inequality into PSD_d.
The face of the PPT cone at W is controlled by rank(W^Gamma), giving rank bounds.

### 2. Lukas: Nonlinear Constraints (Manifold Framework)

**Idea:** If equality constraints involve nonlinear functions and the target value is
a regular value, the constraint set is a manifold. The perturbation direction D lives
in the kernel of the differential of the restricted map (tangent space). The rank/dimension
constraint holds for the differential.

**Open questions:**
- Does the flow in the kernel direction reach a face of lower dimension?
  (This is the analogue of Martin's Rockafellar citation for the linear case.)
- Are there applications of this nonlinear extension?

### 3. Alex: Obstruction to Nonlinear Localisation

**Observation:** The BP argument cannot be localised for nonlinear constraints.

**Counterexample idea:** A nonlinear constraint can force a point to lie deep inside a
convex set (e.g., on a shell strictly inside a polyhedron). In this case, the point is
in the interior of a high-dimensional face, and no perturbation in the kernel will reach
the boundary.

**Key point:** Linear constraints avoid this because they generate suitable faces upon
intersection with the cone. This is exactly the face-generation property used in Step 2
of the proof: the perturbation x + tD must eventually leave V+ (by pointedness), and
x* = x + t+ D lies on a proper face.

### 4. Timo: Inactive Inequality Constraints

**Observation:** If there exists a strictly feasible solution with k inactive constraints,
then d_F <= m - k. Inactive constraints impose no local restriction on the
perturbation direction.

**Proof:** For an inactive constraint alpha_j(x) < b_j, any perturbation D satisfies
alpha_j(x + tD) < b_j for sufficiently small |t|, regardless of alpha_j(D). So
inactive constraints don't contribute to the kernel condition.

### 5. Lennart: Slack Variable Analysis (Equality vs. Inequality)

**Concrete example:** PSD_n with m scalar inequality constraints tr(A_i X) <= b_i.

**Slack formulation:** Introduce s_i >= 0 with tr(A_i X) + s_i = b_i.
This is a problem over the product cone PSD_n x R_+^m with m equality constraints.

**Face dimension in product cone:** dim(face(X, s)) = r(r+1)/2 + |{s_i > 0}|.
BP gives: r(r+1)/2 + |{s_i > 0}| <= m.

**Practical issue:** The iterative face-reduction procedure cannot control whether it
reduces rank(X) or sets slack variables s_i to zero. Without extra structure, we may
just make slacks inactive without reducing the rank of the PSD variable.

**For conic inequality slacks:** If Phi(X) - B = S with S in PSD_k, then
dim(face(S in PSD_k)) = r_S(r_S+1)/2. The slack approach gives:
r_X(r_X+1)/2 + r_S(r_S+1)/2 <= k(k+1)/2 + m_scalar.
This is tighter than Martin's d_F <= k(k+1)/2 + m_scalar when r_S > 0.

### 6. Timo: SOS Relaxations and Polynomial Matrix Inequalities

**Question:** Can the framework handle polynomial matrix inequalities (connecting to
SOS relaxations)? The feasible set is generally no longer convex.

**Status:** Open. Convexity is essential for the face-generation argument (Step 2).
Non-convex sets don't have faces in the usual sense.

## Numerical Verification: Spectrahedra Examples

See `numerics/spectrahedra_bp.jl` for Julia code implementing the following examples.

### Example 1: Block-Diagonal LMI Spectrahedron

Two 4x4 LMI blocks with 10 free variables and 2 trace equality constraints.
- Individual BP: r_F <= 4, r_G <= 4 (weak)
- Block-diagonal BP: r_F + r_G <= 4 (tight)
- Numerical trials confirm BP at extreme points; interior-point solver returns
  analytic centre with higher combined rank.

### Example 2: Lennart's Inequality Constraints (Verified)

PSD_6 with m scalar inequality constraints:
- m=15: rank=2, active=2, inactive=13, face_dim=16 <= 16 (BP tight!)
- Timo's inactive-constraint reduction verified: effective bound uses only active constraints.

### Example 3: PPT Entanglement Witness

Optimising over PPT cone {W >= 0 : W^Gamma >= 0} for 2x2 system:
- m=1 (trace only): rank(W) = rank(W^Gamma) = 2
- m=8: rank(W) = 4 (full), rank(W^Gamma) = 3
- The PPT constraint (conic inequality) bounds rank(W^Gamma), not rank(W).

### Example 4: Compression Constraint

8x8 PSD with 3x3 subblock constraint P'XP >= B:
- Without LMI: BP gives tight scalar bound
- With LMI: Martin's bound correctly accounts for conic constraint dimension
- Lennart's slack approach gives slightly tighter bound when slack is non-zero

## Connections and Open Questions

1. **Martin's spectrahedra BP:** The conic inequality framework naturally handles
   spectrahedra defined by multiple LMIs. The block-diagonal formulation gives joint
   rank bounds strictly tighter than treating each LMI separately.

2. **Commutant interpretation:** For C*-algebra state spaces with conic constraints
   (e.g., PPT), the face-commutant correspondence translates the conic BP bound into
   constraints on the commutant structure. This could characterise which commutant
   decompositions are compatible with PPT-like constraints.

3. **Entanglement witness rank:** Martin's bound gives rank(W^Gamma) <= sqrt(m+1) for
   a PPT witness under m scalar constraints. Combined with the relationship between
   W and W^Gamma, this constrains the Schmidt number of detectable entanglement.

4. **Practical face reduction:** Lennart's observation about uncontrollable slack
   reduction is a genuine computational obstacle. For conic constraints, the slack
   structure (rank of S = Phi(X) - B) determines the tightness of the bound but
   cannot be prescribed a priori.

5. **Nonlinear constraints:** Alex's obstruction is fundamental. The BP argument
   requires the perturbation to reach the boundary of the cone, which fails for
   nonlinear constraints that create "interior shells." Linear constraints are special
   because they generate faces upon intersection.
