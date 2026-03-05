#!/usr/bin/env julia
#
# marginals_bp.jl — Barvinok–Pataki for the quantum marginal problem
#
# Demonstrates:
#   1. Quantum marginal problem as spectrahedron: ρ_AB ≥ 0, Tr_B(ρ) = ρ_A, Tr_A(ρ) = ρ_B
#   2. BP bound: r(r+1)/2 ≤ d_A(d_A+1)/2 + d_B(d_B+1)/2 - 1 for extreme compatible states (real SDP)
#   3. Numerical verification: qubit-qubit, qubit-qutrit, qutrit-qutrit
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

"""Random density matrix of given rank."""
function randdm(n, r)
    V = randn(n, r)
    ρ = V * V'
    return Symmetric(ρ / tr(ρ))
end

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

"""Partial trace over subsystem A for d_A × d_B system."""
function partial_trace_A(X, dA, dB)
    ρB = zeros(eltype(X), dB, dB)
    for ib in 1:dB, jb in 1:dB
        for ia in 1:dA
            ρB[ib, jb] += X[(ia-1)*dB+ib, (ia-1)*dB+jb]
        end
    end
    return ρB
end


# ============================================================================
# Marginal-constrained density matrix SDP
# ============================================================================
#
# Given marginals ρ_A ∈ PSD_{d_A} and ρ_B ∈ PSD_{d_B}, find a joint state
# ρ_AB ∈ PSD_{d_A·d_B} with:
#   Tr_B(ρ_AB) = ρ_A    (d_A(d_A+1)/2 real constraints)
#   Tr_A(ρ_AB) = ρ_B    (d_B(d_B+1)/2 real constraints)
#
# The constraint Tr(ρ_AB) = 1 is already implied by Tr(ρ_A) = 1.
# Total independent constraints:
#   m = d_A(d_A+1)/2 + d_B(d_B+1)/2 - 1
# (subtract 1 because Tr(ρ_A) = Tr(ρ_B) = 1 gives one redundant trace equation)
#
# BP bound: r(r+1)/2 ≤ m  where r = rank(ρ_AB)

"""
Solve the marginal-constrained SDP for given marginals.

Returns (rank, m_total, bp_bound).
"""
function marginal_sdp(dA, dB, ρ_A, ρ_B; n_trials=10)
    d = dA * dB

    # Constraint count
    m_A = dA * (dA + 1) ÷ 2
    m_B = dB * (dB + 1) ÷ 2
    m_total = m_A + m_B - 1  # -1 for redundant trace
    bp_bound = bp_rank_real(m_total)

    Random.seed!(42)
    max_rank = 0

    for trial in 1:n_trials
        model = make_model()
        @variable(model, ρ[1:d, 1:d], PSD)

        # Marginal A constraint: Tr_B(ρ) = ρ_A
        for ia in 1:dA, ja in ia:dA
            lhs = @expression(model,
                sum(ρ[(ia-1)*dB+ib, (ja-1)*dB+ib] for ib in 1:dB))
            @constraint(model, lhs == ρ_A[ia, ja])
        end

        # Marginal B constraint: Tr_A(ρ) = ρ_B
        for ib in 1:dB, jb in ib:dB
            lhs = @expression(model,
                sum(ρ[(ia-1)*dB+ib, (ia-1)*dB+jb] for ia in 1:dA))
            @constraint(model, lhs == ρ_B[ib, jb])
        end

        # Random objective to push to extreme point
        C = randsym(d)
        @objective(model, Min, sum(C[i,j] * ρ[i,j] for i in 1:d, j in 1:d))

        optimize!(model)
        if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
            ρval = value.(ρ)
            ρval = (ρval + ρval') / 2
            r = numrank(ρval)
            max_rank = max(max_rank, r)

            # Verify marginals
            if trial == 1
                ρA_check = partial_trace_B(ρval, dA, dB)
                ρB_check = partial_trace_A(ρval, dA, dB)
                err_A = norm(ρA_check - ρ_A)
                err_B = norm(ρB_check - ρ_B)
                if err_A > 1e-4 || err_B > 1e-4
                    @printf("    Warning: marginal mismatch A=%.2e, B=%.2e\n", err_A, err_B)
                end
            end
        end
    end

    return max_rank, m_total, bp_bound
end


# ============================================================================
# Example 1: Qubit-qubit (2×2)
# ============================================================================
function example1_qubit_qubit()
    println("\n" * "="^72)
    println("Example 1: Qubit-qubit marginal problem (d_A=2, d_B=2)")
    println("="^72)

    dA, dB = 2, 2
    d = dA * dB
    m_A = dA * (dA + 1) ÷ 2  # = 3
    m_B = dB * (dB + 1) ÷ 2  # = 3
    m_total = m_A + m_B - 1   # = 5
    bp_bound = bp_rank_real(m_total)

    println("\n  d_A=$dA, d_B=$dB, d=$d")
    println("  m_A=$m_A (Tr_B constraint), m_B=$m_B (Tr_A constraint)")
    println("  m_total=$m_total (−1 for redundant trace)")
    println("  BP bound: r ≤ $bp_bound")

    Random.seed!(123)

    # Test with various marginals
    test_cases = [
        ("maximally mixed", Matrix{Float64}(I/2, 2, 2), Matrix{Float64}(I/2, 2, 2)),
        ("pure A, mixed B", [1.0 0; 0 0], Matrix{Float64}(I/2, 2, 2)),
        ("both pure", [1.0 0; 0 0], [1.0 0; 0 0]),
        ("random", randdm(2, 2), randdm(2, 2)),
    ]

    println("\n  Marginals             | rank | m   | BP bound | ok")
    println("  " * "-"^60)

    for (name, ρA, ρB) in test_cases
        r, m, bp = marginal_sdp(dA, dB, ρA, ρB; n_trials=10)
        ok = r <= bp
        @printf("  %-22s | %d    | %d   | ≤ %d      | %s\n",
                name, r, m, bp, ok ? "✓" : "✗")
    end
end


# ============================================================================
# Example 2: Qubit-qutrit (2×3)
# ============================================================================
function example2_qubit_qutrit()
    println("\n" * "="^72)
    println("Example 2: Qubit-qutrit marginal problem (d_A=2, d_B=3)")
    println("="^72)

    dA, dB = 2, 3
    d = dA * dB
    m_A = dA * (dA + 1) ÷ 2  # = 3
    m_B = dB * (dB + 1) ÷ 2  # = 6
    m_total = m_A + m_B - 1   # = 8
    bp_bound = bp_rank_real(m_total)

    println("\n  d_A=$dA, d_B=$dB, d=$d")
    println("  m_A=$m_A, m_B=$m_B, m_total=$m_total")
    println("  BP bound: r ≤ $bp_bound")

    Random.seed!(456)

    test_cases = [
        ("maximally mixed", Matrix{Float64}(I/2, 2, 2), Matrix{Float64}(I/3, 3, 3)),
        ("random", randdm(2, 2), randdm(3, 3)),
        ("pure A, mixed B", [1.0 0; 0 0], randdm(3, 2)),
    ]

    println("\n  Marginals             | rank | m   | BP bound | ok")
    println("  " * "-"^60)

    for (name, ρA, ρB) in test_cases
        r, m, bp = marginal_sdp(dA, dB, ρA, ρB; n_trials=10)
        ok = r <= bp
        @printf("  %-22s | %d    | %d   | ≤ %d      | %s\n",
                name, r, m, bp, ok ? "✓" : "✗")
    end
end


# ============================================================================
# Example 3: Qutrit-qutrit (3×3)
# ============================================================================
function example3_qutrit_qutrit()
    println("\n" * "="^72)
    println("Example 3: Qutrit-qutrit marginal problem (d_A=3, d_B=3)")
    println("="^72)

    dA, dB = 3, 3
    d = dA * dB
    m_A = dA * (dA + 1) ÷ 2  # = 6
    m_B = dB * (dB + 1) ÷ 2  # = 6
    m_total = m_A + m_B - 1   # = 11
    bp_bound = bp_rank_real(m_total)

    println("\n  d_A=$dA, d_B=$dB, d=$d")
    println("  m_A=$m_A, m_B=$m_B, m_total=$m_total")
    println("  BP bound: r ≤ $bp_bound")

    Random.seed!(789)

    test_cases = [
        ("maximally mixed", Matrix{Float64}(I/3, 3, 3), Matrix{Float64}(I/3, 3, 3)),
        ("random", randdm(3, 3), randdm(3, 3)),
        ("random rank-2", randdm(3, 2), randdm(3, 2)),
    ]

    println("\n  Marginals             | rank | m   | BP bound | ok")
    println("  " * "-"^60)

    for (name, ρA, ρB) in test_cases
        r, m, bp = marginal_sdp(dA, dB, ρA, ρB; n_trials=10)
        ok = r <= bp
        @printf("  %-22s | %d    | %d   | ≤ %d      | %s\n",
                name, r, m, bp, ok ? "✓" : "✗")
    end
end


# ============================================================================
# Summary table
# ============================================================================
function summary_table()
    println("\n" * "="^72)
    println("Summary: BP bounds for quantum marginal problem")
    println("="^72)

    println("\n  d_A  d_B  d=d_A·d_B  m_A      m_B      m_total  BP: r ≤")
    println("  " * "-"^60)

    for (dA, dB) in [(2,2), (2,3), (2,4), (3,3), (3,4), (4,4)]
        d = dA * dB
        m_A = dA * (dA + 1) ÷ 2
        m_B = dB * (dB + 1) ÷ 2
        m_total = m_A + m_B - 1
        bp_bound = bp_rank_real(m_total)
        @printf("  %d    %d    %2d         %2d       %2d       %2d       %d\n",
                dA, dB, d, m_A, m_B, m_total, bp_bound)
    end
end


# ============================================================================
# Run all examples
# ============================================================================

println("Barvinok-Pataki for the Quantum Marginal Problem")
println("=" ^ 50)
println("Solver: $solver_available")

example1_qubit_qubit()
example2_qubit_qutrit()
example3_qutrit_qutrit()
summary_table()

println("\n" * "="^72)
println("Summary:")
println("="^72)
println("""
1. The quantum marginal problem (finding ρ_AB compatible with given ρ_A, ρ_B)
   is a spectrahedron: ρ_AB ≥ 0 with linear marginal constraints.

2. Constraint count: m = d_A(d_A+1)/2 + d_B(d_B+1)/2 - 1
   (real symmetric SDP, subtract 1 for redundant trace normalisation).

3. BP bound: r(r+1)/2 ≤ m, i.e., rank(ρ_AB) ≤ floor((-1+√(1+8m))/2).
   For qubits (2×2): m=5, r ≤ 2.
   For qubit-qutrit (2×3): m=8, r ≤ 3.
   For qutrits (3×3): m=11, r ≤ 4.

4. This means extreme compatible states always have low rank relative to
   the ambient dimension d_A·d_B, controlling the entanglement structure.
""")
