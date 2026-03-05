#!/usr/bin/env julia
#
# thermal_bp.jl — Barvinok–Pataki for thermal operations (Gibbs-preserving channels)
#
# Demonstrates:
#   1. Thermal operations as spectrahedron: CPTP + Φ(γ_β) = γ_β
#   2. BP bound on Kraus complexity of extreme thermal operations
#   3. Temperature dependence: β→0 recovers doubly stochastic, β→∞ is ground-state preserving
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

"""Random real symmetric matrix."""
randsym(n) = let M = randn(n, n); (M + M') / 2 end

"""Partial trace over subsystem B for d_A × d_B system."""
function partial_trace_B(X, dA, dB)
    ρA = zeros(eltype(X), dA, dA)
    for ia in 1:dA, ja in 1:dA
        for ib in 1:dB
            ρA[ia, ja] += X[(ia-1)*dB+ib, (ja-1)*dB+ib]
        end
    end
    return ρA
end

"""
Gibbs state γ_β = exp(-β H) / Z for Hamiltonian H with eigenvalues E.

Returns a diagonal density matrix.
"""
function gibbs_state(E, β)
    d = length(E)
    if β == 0
        return Matrix{Float64}(I/d, d, d)  # maximally mixed
    end
    w = exp.(-β .* E)
    w ./= sum(w)
    return diagm(w)
end


# ============================================================================
# Gibbs-preserving channel SDP
# ============================================================================
#
# A Gibbs-preserving channel Φ: M_d → M_d satisfies:
#   1. CPTP: J_Φ ∈ PSD_{d²}, Tr_out(J_Φ) = I_d
#   2. Gibbs preservation: Φ(γ_β) = γ_β
#
# In Choi picture, Φ(γ_β) = γ_β translates to:
#   Tr_in((γ_β^T ⊗ I_out) J_Φ) = γ_β
# which gives d(d+1)/2 real constraints (upper triangle of d×d Hermitian).
#
# Total constraints: m = m_TP + m_Gibbs
#   m_TP = d(d+1)/2
#   m_Gibbs = d(d+1)/2
# But some constraints may overlap (e.g., trace(Φ(γ_β)) = 1 is implied by TP).
# Conservatively: m = d(d+1)/2 + d(d+1)/2 - 1 = d(d+1) - 1

"""
Build and solve thermal operation SDP for given Hamiltonian eigenvalues and β.

Returns (max_rank, m_total, bp_bound).
"""
function thermal_sdp(E, β; n_trials=10)
    d = length(E)
    d_choi = d^2
    γ = gibbs_state(E, β)

    # Constraint counts
    m_tp = d * (d + 1) ÷ 2
    m_gibbs = d * (d + 1) ÷ 2
    m_total = m_tp + m_gibbs - 1  # -1 for redundant trace
    bp_bound = bp_rank_real(m_total)

    Random.seed!(42)
    max_rank = 0

    for trial in 1:n_trials
        model = make_model()
        @variable(model, J[1:d_choi, 1:d_choi], PSD)

        # TP constraint: Tr_out(J) = I_d
        for ia in 1:d, ja in ia:d
            lhs = @expression(model,
                sum(J[(ia-1)*d + ib, (ja-1)*d + ib] for ib in 1:d))
            rhs = (ia == ja) ? 1.0 : 0.0
            @constraint(model, lhs == rhs)
        end

        # Gibbs preservation: Tr_in((γ^T ⊗ I) J) = γ
        # [Tr_in((γ^T ⊗ I) J)]_{ib,jb} = Σ_{ia} γ[ia,ia] J[(ia-1)*d+ib, (ia-1)*d+jb]
        # (γ is diagonal, so γ^T = γ)
        for ib in 1:d, jb in ib:d
            lhs = @expression(model,
                sum(γ[ia, ia] * J[(ia-1)*d + ib, (ia-1)*d + jb] for ia in 1:d))
            rhs = γ[ib, jb]
            @constraint(model, lhs == rhs)
        end

        # Random objective
        C = randsym(d_choi)
        @objective(model, Min, sum(C[i,j] * J[i,j] for i in 1:d_choi, j in 1:d_choi))

        optimize!(model)
        if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
            Jval = value.(J)
            Jval = (Jval + Jval') / 2
            r = numrank(Jval)
            max_rank = max(max_rank, r)
        end
    end

    return max_rank, m_total, bp_bound
end


# ============================================================================
# Example 1: Qubit thermal operations at various temperatures
# ============================================================================
function example1_qubit_thermal()
    println("\n" * "="^72)
    println("Example 1: Qubit Gibbs-preserving channels (d=2)")
    println("="^72)

    d = 2
    E = [0.0, 1.0]  # qubit Hamiltonian with gap Δ=1

    m_tp = d * (d + 1) ÷ 2  # = 3
    m_gibbs = d * (d + 1) ÷ 2  # = 3
    m_total = m_tp + m_gibbs - 1  # = 5
    bp_bound = bp_rank_real(m_total)

    println("\n  d=$d, m_TP=$m_tp, m_Gibbs=$m_gibbs, m_total=$m_total")
    println("  BP bound: r ≤ $bp_bound")

    βs = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    println("\n  β       γ_ground    γ_excited   rank   BP bound   ok")
    println("  " * "-"^60)

    for β in βs
        γ = gibbs_state(E, β)
        r, m, bp = thermal_sdp(E, β; n_trials=10)
        ok = r <= bp
        @printf("  %-7.1f  %.4f      %.4f       %d      ≤ %d        %s\n",
                β, γ[1,1], γ[2,2], r, bp, ok ? "✓" : "✗")
    end
end


# ============================================================================
# Example 2: Qutrit thermal operations
# ============================================================================
function example2_qutrit_thermal()
    println("\n" * "="^72)
    println("Example 2: Qutrit Gibbs-preserving channels (d=3)")
    println("="^72)

    d = 3
    E = [0.0, 1.0, 2.5]  # qutrit with non-degenerate spectrum

    m_tp = d * (d + 1) ÷ 2  # = 6
    m_gibbs = d * (d + 1) ÷ 2  # = 6
    m_total = m_tp + m_gibbs - 1  # = 11
    bp_bound = bp_rank_real(m_total)

    println("\n  d=$d, E=$E")
    println("  m_TP=$m_tp, m_Gibbs=$m_gibbs, m_total=$m_total")
    println("  BP bound: r ≤ $bp_bound")

    βs = [0.0, 0.5, 1.0, 2.0, 5.0]

    println("\n  β       γ₁      γ₂      γ₃      rank   BP bound   ok")
    println("  " * "-"^65)

    for β in βs
        γ = gibbs_state(E, β)
        r, m, bp = thermal_sdp(E, β; n_trials=10)
        ok = r <= bp
        @printf("  %-7.1f  %.3f   %.3f   %.3f    %d      ≤ %d        %s\n",
                β, γ[1,1], γ[2,2], γ[3,3], r, bp, ok ? "✓" : "✗")
    end
end


# ============================================================================
# Example 3: Temperature dependence — qubit BP bound vs β
# ============================================================================
function example3_temperature_dependence()
    println("\n" * "="^72)
    println("Example 3: Temperature dependence of BP bound (qubit)")
    println("="^72)

    d = 2
    E = [0.0, 1.0]

    # At β=0: γ = I/d (maximally mixed), Gibbs preservation = doubly stochastic
    # At β→∞: γ → |0⟩⟨0|, ground-state preservation adds a rank-1 constraint
    # In both limits, the constraint count is the same (m = d(d+1) - 1),
    # but the effective constraints change character.

    println("\n  β=0 limit: Gibbs preservation = unital (doubly stochastic)")
    println("  β→∞ limit: ground-state preserving channel")
    println()

    # The constraint count m doesn't depend on β (same number of equations),
    # but the tightness of the BP bound may vary because some constraints
    # become linearly dependent at special temperatures.

    βs = [0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]

    println("  β         m_total  BP bound  observed rank  gap")
    println("  " * "-"^55)

    for β in βs
        r, m, bp = thermal_sdp(E, β; n_trials=10)
        @printf("  %-9.2f  %d        ≤ %d       %d              %d\n",
                β, m, bp, r, bp - r)
    end
end


# ============================================================================
# Summary table
# ============================================================================
function summary_table()
    println("\n" * "="^72)
    println("Summary: BP bounds for Gibbs-preserving channels")
    println("="^72)

    println("\n  d   m_TP     m_Gibbs  m_total  BP: r ≤   (at β=1)")
    println("  " * "-"^55)

    for d in [2, 3, 4, 5]
        E = collect(0.0:1.0:(d-1))
        m_tp = d * (d + 1) ÷ 2
        m_gibbs = d * (d + 1) ÷ 2
        m_total = m_tp + m_gibbs - 1
        bp_bound = bp_rank_real(m_total)

        r, _, _ = thermal_sdp(E, 1.0; n_trials=5)

        @printf("  %d   %3d      %3d      %3d      %3d        r_obs=%d\n",
                d, m_tp, m_gibbs, m_total, bp_bound, r)
    end
end


# ============================================================================
# Run all examples
# ============================================================================

println("Barvinok-Pataki for Thermal Operations")
println("=" ^ 50)
println("Solver: $solver_available")

example1_qubit_thermal()
example2_qutrit_thermal()
example3_temperature_dependence()
summary_table()

println("\n" * "="^72)
println("Summary:")
println("="^72)
println("""
1. Thermal operations (Gibbs-preserving CPTP maps) form a spectrahedron:
   J_Φ ≥ 0, Tr_out(J) = I (TP), and Φ(γ_β) = γ_β (Gibbs preservation).

2. Constraint count: m = d(d+1)/2 (TP) + d(d+1)/2 (Gibbs) - 1 = d(d+1) - 1.
   BP bound: r(r+1)/2 ≤ m.

3. Temperature dependence:
   - β → 0: γ_β → I/d, Gibbs preservation becomes unitality (doubly stochastic).
   - β → ∞: γ_β → |0⟩⟨0|, Gibbs preservation becomes ground-state fixing.
   The constraint count is the same, but the effective rank of extremals may vary.

4. The BP bound gives the maximum Kraus complexity of extreme thermal operations,
   relevant for resource theories of thermodynamics.
""")
