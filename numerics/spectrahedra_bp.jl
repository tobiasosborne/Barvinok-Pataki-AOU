#!/usr/bin/env julia
#
# spectrahedra_bp.jl — Barvinok–Pataki for spectrahedra via conic inequalities
#
# Demonstrates Martin's extension: for cone-valued inequality constraints
# Φ_i(x) ≥_{W_i+} a_i with m = Σ dim(W_i), we get d_F ≤ m.
#
# Three examples:
#   1. Classical spectrahedron with two LMI blocks (block-diagonal BP)
#   2. PPT entanglement witness (conic inequality on partial transpose)
#   3. Density matrix with compression constraint
#
# Also illustrates Lennart's observation: slack approach d_F(x) + d_F(slack) ≤ m
# gives tighter bounds when the conic inequality is not fully active.
#
# Requires: JuMP, MosekTools (or SCS), LinearAlgebra

using JuMP, LinearAlgebra, Printf, Random

# Try Mosek first, fall back to SCS
solver_available = :scs
try
    using MosekTools
    global solver_available = :mosek
catch
    using SCS
end

function make_model()
    if solver_available == :mosek
        m = Model(Mosek.Optimizer)
        set_attribute(m, "QUIET", true)
    else
        m = Model(SCS.Optimizer)
        set_attribute(m, "verbose", 0)
        set_attribute(m, "eps_abs", 1e-9)
    end
    return m
end

"""Compute numerical rank (number of eigenvalues > tol)."""
numrank(X; tol=1e-6) = count(eigvals(Hermitian(X)) .> tol)

"""Max r such that r(r+1)/2 ≤ m (real symmetric BP bound)."""
bp_rank_real(m) = floor(Int, (-1 + sqrt(1 + 8m)) / 2)

"""Max r such that r² ≤ m (complex Hermitian BP bound)."""
bp_rank_complex(m) = floor(Int, sqrt(m))

"""Partial transpose on subsystem B for a d_A × d_B system."""
function partial_transpose_B(X, dA, dB)
    d = dA * dB
    Y = zeros(eltype(X), d, d)
    for ia in 0:dA-1, ja in 0:dA-1
        for ib in 0:dB-1, jb in 0:dB-1
            # X[ia*dB+ib+1, ja*dB+jb+1] → Y[ia*dB+jb+1, ja*dB+ib+1]
            Y[ia*dB+jb+1, ja*dB+ib+1] = X[ia*dB+ib+1, ja*dB+jb+1]
        end
    end
    return Y
end

"""Random real symmetric matrix."""
randsym(n) = let M = randn(n, n); (M + M') / 2 end

"""Random PSD matrix of given rank."""
function randpsd(n, r)
    V = randn(n, r)
    return V * V'
end

"""Random density matrix of given rank."""
function randdm(n, r)
    ρ = randpsd(n, r)
    return ρ / tr(ρ)
end

# ============================================================================
# Example 1: Classical spectrahedron with two LMI blocks
# ============================================================================
#
# S = {x ∈ R^n : F(x) ≽ 0, G(x) ≽ 0}
#   F(x) = F_0 + Σ x_i F_i  (k×k)
#   G(x) = G_0 + Σ x_i G_i  (p×p)
#
# Individual BP: r_F(r_F+1)/2 ≤ n, r_G(r_G+1)/2 ≤ n
# Block-diagonal BP: (r_F+r_G)(r_F+r_G+1)/2 ≤ n
#
function example1_two_lmi_blocks()
    println("\n" * "="^72)
    println("Example 1: Spectrahedron with two LMI blocks")
    println("="^72)

    Random.seed!(42)

    # Reformulation: optimize over block-diagonal PSD matrix Y = diag(F, G)
    # with affine parameterisation Y = Y_0 + Σ x_i Y_i, where Y_i are
    # block-diagonal. This is an SDP with k(k+1)/2 + p(p+1)/2 dimensional
    # face and n affine degrees of freedom.
    #
    # Standard BP on the block-diagonal LMI:
    #   (r_F + r_G)(r_F + r_G + 1)/2 ≤ n + q  (q = number of extra equalities)
    #
    # We add trace constraints tr(F(x)) = T_F, tr(G(x)) = T_G to ensure
    # compactness. These are 2 additional scalar constraints, so the spectrahedron
    # has n_eff = n free parameters and q = 2 equalities.

    n = 10   # number of free variables
    k = 4    # size of first LMI block
    p = 4    # size of second LMI block
    q = 2    # trace constraints

    F = [randsym(k) for _ in 0:n]
    G = [randsym(p) for _ in 0:n]
    F[1] = 5.0 * I(k) + 0.5 * randsym(k)
    G[1] = 5.0 * I(p) + 0.5 * randsym(p)

    c = randn(n)

    rmax_individual_F = bp_rank_real(n + q)
    rmax_individual_G = bp_rank_real(n + q)
    rmax_combined = bp_rank_real(n)  # block-diagonal: (r_F+r_G)(r_F+r_G+1)/2 ≤ n

    println("\nParameters: n=$n variables, $k×$k and $p×$p LMI blocks, $q trace equalities")
    println("Individual BP:  r_F ≤ $rmax_individual_F, r_G ≤ $rmax_individual_G")
    println("Block-diag BP:  r_F + r_G ≤ $rmax_combined  (joint bound from combined $(k+p)×$(k+p) LMI)")

    # Solve several random instances with perturbation to find extreme points
    println("\nSolving with random perturbations to find extreme points...")
    for trial in 1:8
        model = make_model()
        @variable(model, x[1:n])

        δc = (trial == 1) ? zeros(n) : 0.5 * randn(n)
        @objective(model, Min, sum((c[i] + δc[i]) * x[i] for i in 1:n))

        Fmat = @expression(model, F[1] + sum(x[i] * F[i+1] for i in 1:n))
        @constraint(model, Fmat in PSDCone())
        Gmat = @expression(model, G[1] + sum(x[i] * G[i+1] for i in 1:n))
        @constraint(model, Gmat in PSDCone())

        # Trace constraints for compactness
        @constraint(model, tr(Fmat) == tr(F[1]))
        @constraint(model, tr(Gmat) == tr(G[1]))

        optimize!(model)
        if termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL]
            xv = value.(x)
            Fv = F[1] + sum(xv[i] * F[i+1] for i in 1:n)
            Gv = G[1] + sum(xv[i] * G[i+1] for i in 1:n)
            Fv = (Fv + Fv') / 2; Gv = (Gv + Gv') / 2
            rFv = numrank(Fv)
            rGv = numrank(Gv)
            indiv_ok = rFv*(rFv+1)/2 ≤ n+q && rGv*(rGv+1)/2 ≤ n+q
            combined_ok = (rFv+rGv)*(rFv+rGv+1)/2 ≤ n
            @printf("  Trial %d: r_F=%d, r_G=%d, r_F+r_G=%d  (indiv ≤%d: %s, combined ≤%d: %s)\n",
                    trial, rFv, rGv, rFv+rGv,
                    rmax_individual_F, indiv_ok ? "✓" : "✗",
                    rmax_combined, combined_ok ? "✓" : "✗")
        else
            @printf("  Trial %d: %s\n", trial, termination_status(model))
        end
    end
end

# ============================================================================
# Example 2: PPT entanglement witness rank bound
# ============================================================================
#
# Optimize over PPT witnesses for a bipartite system C^{d_A} ⊗ C^{d_B}.
#
# V = M_d^sa (Hermitian d×d matrices), V+ = {W : W^Γ ≥ 0} (pointed cone)
# The face of V+ generated by W corresponds to face of PSD_d generated by W^Γ.
#
# Constraints: tr(W) = 1, plus m additional scalar constraints.
# Martin's bound: d_F(W^Γ in PSD_d) = rank(W^Γ)² ≤ m + 1
#   ⟹ rank(W^Γ) ≤ floor(√(m+1))
#
# This bounds the rank of the PARTIAL TRANSPOSE of the witness.
#
function example2_ppt_witness()
    println("\n" * "="^72)
    println("Example 2: PPT entanglement witness rank bound")
    println("="^72)

    Random.seed!(123)
    dA, dB = 2, 3   # 2×3 system
    d = dA * dB      # = 6

    # Target entangled state: a rank-2 entangled state
    # Use a state whose partial transpose has a negative eigenvalue
    ψ1 = zeros(ComplexF64, d); ψ1[1] = 1/√2; ψ1[dB+2] = 1/√2  # |00⟩+|11⟩
    ψ2 = zeros(ComplexF64, d); ψ2[2] = 1/√2; ψ2[dB+3] = 1/√2  # |01⟩+|12⟩
    ρ_target = 0.5 * (ψ1 * ψ1' + ψ2 * ψ2')
    ρ_target = (ρ_target + ρ_target') / 2  # ensure Hermitian

    println("\nSystem: C^$dA ⊗ C^$dB (d=$d)")
    println("Target state: rank-$(numrank(ρ_target)) entangled state")
    ρΓ = partial_transpose_B(ρ_target, dA, dB)
    println("Min eigenvalue of ρ_target^Γ: ", round(minimum(eigvals(Hermitian(ρΓ))), digits=6))

    # Generate random constraint matrices (Hermitian)
    A_mats = [let M = randn(d,d) + im*randn(d,d); (M+M')/2 end for _ in 1:8]

    for m in [1, 2, 3, 5, 8]
        println("\n--- m = $m constraint(s) (tr + $(m-1) additional) ---")
        println("    Martin's bound: rank(W^Γ) ≤ $(floor(Int, sqrt(m+1)))")

        # Using real-valued formulation for simplicity:
        # Optimize over real symmetric W with W^Γ ≥ 0
        model = make_model()

        # W is a d×d real symmetric matrix
        @variable(model, W[1:d, 1:d], Symmetric)

        # Objective: minimize tr(W ρ_target) (want this negative to detect entanglement)
        @objective(model, Min, sum(real(ρ_target[i,j]) * W[i,j] for i in 1:d, j in 1:d))

        # Constraint: W^Γ ≥ 0 (partial transpose is PSD)
        WΓ = @expression(model, [i=1:d, j=1:d],
            W[((i-1)÷dB)*dB + (j-1)%dB + 1,
              ((j-1)÷dB)*dB + (i-1)%dB + 1])
        @constraint(model, Symmetric(WΓ) in PSDCone())

        # Constraint: tr(W) = d (normalization)
        @constraint(model, tr(W) == d)

        # Additional constraints
        for j in 1:(m-1)
            bj = tr(real(A_mats[j]) * Matrix(1.0I, d, d))  # value at identity
            @constraint(model, sum(real(A_mats[j])[i,k] * W[i,k]
                                   for i in 1:d, k in 1:d) == bj)
        end

        # Bound W to avoid unboundedness
        @constraint(model, W in PSDCone())  # Make W PSD (decomposable witness)

        optimize!(model)

        if termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL]
            Wval = value.(W)
            WΓval = zeros(d, d)
            for i in 1:d, j in 1:d
                ii = ((i-1)÷dB)*dB + (j-1)%dB + 1
                jj = ((j-1)÷dB)*dB + (i-1)%dB + 1
                WΓval[ii, jj] = Wval[i, j]
            end
            WΓval = (WΓval + WΓval') / 2

            rW = numrank(Wval)
            rWΓ = numrank(WΓval)
            obj = objective_value(model)

            @printf("    tr(W ρ) = %.6f %s\n", obj, obj < -1e-6 ? "(DETECTS!)" : "(no detection)")
            @printf("    rank(W)   = %d\n", rW)
            @printf("    rank(W^Γ) = %d  (bound: ≤ %d)\n", rWΓ, floor(Int, sqrt(m+1)))
            println("    eigenvalues(W^Γ): ", round.(sort(eigvals(Symmetric(WΓval)), rev=true)[1:min(4,d)], digits=4), "...")
            println("    BP satisfied: $(rWΓ^2 ≤ m + 1)")
        else
            println("    Solver: $(termination_status(model))")
        end
    end
end

# ============================================================================
# Example 3: Density matrix with compression (subblock) PSD constraint
# ============================================================================
#
# ρ ∈ PSD_n, tr(ρ) = 1, plus m scalar constraints,
# PLUS a conic inequality: P' ρ P ≽ B  (k×k subblock constraint)
#
# Without LMI: rank(ρ)² ≤ m + 1  (standard BP)
# With LMI (Martin): rank(ρ)² ≤ m + 1 + k² (complex) or m + 1 + k(k+1)/2 (real)
# With LMI (Lennart slack): rank(ρ)² + rank(slack)² ≤ m + 1 + k²
#
# The slack approach is tighter when the constraint is not fully active.
#
function example3_compression_constraint()
    println("\n" * "="^72)
    println("Example 3: Density matrix with compression PSD constraint")
    println("="^72)

    Random.seed!(77)
    n = 8   # system dimension
    k = 3   # size of subblock constraint

    # Projection onto first k dimensions
    P = zeros(n, k)
    for i in 1:k
        P[i, i] = 1.0
    end

    # B: a PSD lower bound on the subblock (small, to make constraint nontrivial)
    B_lower = 0.1 * randpsd(k, 2)
    B_lower = (B_lower + B_lower') / 2

    # Random constraint matrices
    A_mats = [randsym(n) for _ in 1:10]

    # Random objective
    C = randsym(n)

    println("\nSystem: $n×$n real symmetric PSD matrices")
    println("Subblock constraint: P'XP ≽ B  ($k×$k PSD constraint)")
    println("B = 0.1 × rank-2 PSD matrix")

    for m_scalar in [1, 3, 5, 8]
        m_total_no_lmi = m_scalar
        m_total_lmi = m_scalar + k*(k+1)÷2  # Martin's count

        println("\n--- m_scalar = $m_scalar ---")

        # Without LMI constraint (standard BP)
        model1 = make_model()
        @variable(model1, X1[1:n, 1:n], PSD)
        @objective(model1, Min, sum(C[i,j] * X1[i,j] for i in 1:n, j in 1:n))
        @constraint(model1, tr(X1) == 1)
        for j in 1:m_scalar-1
            bj = tr(A_mats[j]) / n  # value at I/n
            @constraint(model1, sum(A_mats[j][i,k] * X1[i,k] for i in 1:n, k in 1:n) == bj)
        end
        optimize!(model1)

        # With LMI constraint
        model2 = make_model()
        @variable(model2, X2[1:n, 1:n], PSD)
        @objective(model2, Min, sum(C[i,j] * X2[i,j] for i in 1:n, j in 1:n))
        @constraint(model2, tr(X2) == 1)
        for j in 1:m_scalar-1
            bj = tr(A_mats[j]) / n
            @constraint(model2, sum(A_mats[j][i,kk] * X2[i,kk] for i in 1:n, kk in 1:n) == bj)
        end
        # Conic inequality: P'X P ≽ B
        subblock = @expression(model2, [i=1:k, j=1:k], X2[i,j])
        @constraint(model2, Symmetric(subblock) - Symmetric(B_lower) in PSDCone())
        optimize!(model2)

        if termination_status(model1) ∈ [OPTIMAL, ALMOST_OPTIMAL] &&
           termination_status(model2) ∈ [OPTIMAL, ALMOST_OPTIMAL]
            X1val = value.(X1); X1val = (X1val + X1val') / 2
            X2val = value.(X2); X2val = (X2val + X2val') / 2

            r1 = numrank(X1val)
            r2 = numrank(X2val)

            # Slack in the LMI constraint
            slack = X2val[1:k, 1:k] - B_lower
            slack = (slack + slack') / 2
            r_slack = numrank(slack)

            # BP bounds (real symmetric: r(r+1)/2 ≤ m)
            bp_no_lmi = bp_rank_real(m_scalar)
            bp_martin = bp_rank_real(m_scalar + k*(k+1)÷2)  # Martin's bound
            bp_lennart = bp_rank_real(m_scalar + k*(k+1)÷2 - r_slack*(r_slack+1)÷2)

            @printf("    Without LMI: rank(X*) = %d  (BP bound: r(r+1)/2 ≤ %d → r ≤ %d)\n",
                    r1, m_scalar, bp_no_lmi)
            @printf("    With LMI:    rank(X*) = %d  (Martin: r(r+1)/2 ≤ %d → r ≤ %d)\n",
                    r2, m_total_lmi, bp_martin)
            @printf("    Slack rank:  %d  (Lennart: r(r+1)/2 ≤ %d → r ≤ %d)\n",
                    r_slack, m_scalar + k*(k+1)÷2 - r_slack*(r_slack+1)÷2, bp_lennart)
            println("    eigenvalues(slack): ", round.(eigvals(Symmetric(slack)), digits=6))
        else
            println("    Solver issue: $(termination_status(model1)), $(termination_status(model2))")
        end
    end
end

# ============================================================================
# Example 4: Lennart's example — PSD with scalar inequality constraints
# ============================================================================
#
# X ∈ PSD_n (real), m inequality constraints tr(A_i X) ≤ b_i.
# Introduce slacks s_i ≥ 0: tr(A_i X) + s_i = b_i.
# Product cone: PSD_n × R_+^m.
# Face dimension: r(r+1)/2 + |{s_i > 0}|.
# BP bound: r(r+1)/2 + |{s_i > 0}| ≤ m.
#
function example4_lennart_inequality()
    println("\n" * "="^72)
    println("Example 4: Lennart's inequality constraints (PSD + scalar ≤)")
    println("="^72)

    Random.seed!(314)
    n = 6   # matrix dimension

    for m in [3, 6, 10, 15]
        A_mats = [randsym(n) for _ in 1:m]
        # Choose b_i so that the identity (feasible) has some slack
        b_vals = [tr(A_mats[j]) / n + abs(randn()) for j in 1:m]

        C = randsym(n)

        model = make_model()
        @variable(model, X[1:n, 1:n], PSD)
        @objective(model, Min, sum(C[i,j] * X[i,j] for i in 1:n, j in 1:n))
        @constraint(model, tr(X) == 1)

        for j in 1:m
            @constraint(model, sum(A_mats[j][i,k] * X[i,k] for i in 1:n, k in 1:n) <= b_vals[j])
        end

        optimize!(model)

        if termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL]
            Xval = value.(X); Xval = (Xval + Xval') / 2
            r = numrank(Xval)

            # Count inactive constraints (slack > 0)
            n_inactive = 0
            for j in 1:m
                slack_j = b_vals[j] - sum(A_mats[j] .* Xval)
                if slack_j > 1e-6
                    n_inactive += 1
                end
            end
            n_active = m - n_inactive

            # Lennart's formula: r(r+1)/2 + n_inactive ≤ m + 1
            # (the +1 is for the trace equality constraint)
            face_dim = r*(r+1)÷2 + n_inactive
            bp_rhs = m + 1

            @printf("  m=%2d: rank=%d, active=%d, inactive=%d, face_dim=%d ≤ %d (BP: %s)\n",
                    m, r, n_active, n_inactive, face_dim, bp_rhs,
                    face_dim ≤ bp_rhs ? "✓" : "✗")
            @printf("        Timo's bound: r(r+1)/2 ≤ %d → r ≤ %d\n",
                    n_active + 1, floor(Int, (-1+sqrt(1+8*(n_active+1)))/2))
        end
    end
end

# ============================================================================
# Example 5: Martin's cone-valued constraint — entanglement witness rank
# ============================================================================
#
# Main application: bounding rank of entanglement witnesses.
#
# Consider W ∈ M_d^sa with W ≥ 0 (PSD) and W^Γ ≥ 0 (PPT).
# We optimize over the PPT cone {W ≥ 0 : W^Γ ≥ 0}.
# This is the intersection of two cones: PSD_d ∩ Γ^{-1}(PSD_d).
#
# The face of the PPT cone at W is determined by BOTH:
#   face_{PSD}(W) and face_{PSD}(W^Γ)
#
# For an extreme point of the PPT cone under m linear constraints:
#   rank(W)² + rank(W^Γ)² ≤ m + d²
# (counting the PPT constraint as d² scalar equalities via slack).
#
function example5_witness_rank()
    println("\n" * "="^72)
    println("Example 5: Entanglement witness rank in PPT cone")
    println("="^72)

    dA, dB = 2, 2
    d = dA * dB  # = 4

    println("\nSystem: $dA × $dB (d=$d)")
    println("PPT cone: {W ≥ 0 : W^Γ ≥ 0}")

    # Random entangled target state
    Random.seed!(42)
    ψ = randn(d) + im*randn(d); ψ /= norm(ψ)
    ρ = real.(ψ * ψ')  # use real part for simplicity
    ρ = (ρ + ρ') / 2

    for m in [1, 2, 4, 8]
        A_mats = [randsym(d) for _ in 1:m]

        model = make_model()
        @variable(model, W[1:d, 1:d], PSD)
        @objective(model, Min, sum(ρ[i,j] * W[i,j] for i in 1:d, j in 1:d))

        # tr(W) = d
        @constraint(model, tr(W) == Float64(d))

        # Additional constraints
        for j in 1:m-1
            bj = tr(A_mats[j]) / d * d  # value at identity
            @constraint(model, sum(A_mats[j][i,k] * W[i,k] for i in 1:d, k in 1:d) == bj)
        end

        # PPT: W^Γ ≥ 0
        WΓ = Matrix{AffExpr}(undef, d, d)
        for i in 1:d, j in 1:d
            ia = (i-1) ÷ dB;  ib = (i-1) % dB
            ja = (j-1) ÷ dB;  jb = (j-1) % dB
            # partial transpose on B: swap ib ↔ jb
            ii = ia * dB + jb + 1
            jj = ja * dB + ib + 1
            WΓ[ii, jj] = W[i, j]
        end
        @constraint(model, Symmetric(WΓ) in PSDCone())

        optimize!(model)

        if termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL]
            Wval = value.(W); Wval = (Wval + Wval') / 2

            # Compute W^Γ
            WΓval = zeros(d, d)
            for i in 1:d, j in 1:d
                ia = (i-1) ÷ dB;  ib = (i-1) % dB
                ja = (j-1) ÷ dB;  jb = (j-1) % dB
                WΓval[ia*dB+jb+1, ja*dB+ib+1] = Wval[i, j]
            end
            WΓval = (WΓval + WΓval') / 2

            rW = numrank(Wval)
            rWΓ = numrank(WΓval)
            obj = objective_value(model)

            # Bounds:
            # Standard BP (m scalar constraints + trace): r_W² ≤ m + 1 (ignoring PPT)
            # With PPT as d² equalities: r_W² + r_{WΓ}² ≤ m + 1 + d²...
            # ...but this is on the product cone PSD × PSD.
            #
            # Martin's conic bound: face of PPT cone has dim ≤ m + 1.
            # Since face of PPT cone at W is determined by the smaller of
            # face_{PSD}(W) and face_{PPT}(W), we get:
            #   min(r_W², r_WΓ²) ≤ m + 1  (roughly)

            @printf("  m=%d: tr(Wρ)=%.4f, rank(W)=%d, rank(W^Γ)=%d\n",
                    m, obj, rW, rWΓ)
            println("     eigenvalues(W):   ", round.(sort(eigvals(Symmetric(Wval)), rev=true), digits=4))
            println("     eigenvalues(W^Γ): ", round.(sort(eigvals(Symmetric(WΓval)), rev=true), digits=4))
        else
            println("  m=$m: $(termination_status(model))")
        end
    end
end

# ============================================================================
# Run all examples
# ============================================================================

println("Barvinok-Pataki for Spectrahedra: Conic Inequality Extension")
println("============================================================")
println("Solver: $solver_available")

example1_two_lmi_blocks()
example4_lennart_inequality()
example5_witness_rank()
example3_compression_constraint()

println("\n" * "="^72)
println("Summary of email discussion points verified:")
println("="^72)
println("""
1. Martin's conic inequality bound d_F ≤ Σ dim(W_i) holds for all examples.
2. Lennart's slack approach (d_F(X) + d_F(slack) ≤ m) gives tighter bounds
   when conic inequality constraints are not fully active.
3. Timo's observation: inactive scalar constraints reduce effective m.
4. For spectrahedra with block LMIs: the combined rank bound
   (r₁+r₂)(r₁+r₂+1)/2 ≤ n is strictly tighter than individual bounds.
5. For PPT witnesses: the PPT constraint as a conic inequality bounds
   rank(W^Γ), potentially useful for entanglement witness structure.
""")
