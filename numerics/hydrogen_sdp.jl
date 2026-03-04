#!/usr/bin/env julia
#
# hydrogen_sdp.jl — Barvinok–Pataki rank bounds for hydrogen in a magnetic field
#
# Solves the SDP:  min Tr(H ρ)  s.t. ρ ≥ 0, Tr(ρ) = 1, Tr(Aⱼ ρ) = bⱼ
# for the hydrogen atom restricted to an energy shell, with angular momentum
# constraints.  Verifies that extreme-point optimisers have rank matching
# the BP bound  r² ≤ m + 1.
#
# Two-phase approach:
#   Phase 1: Solve the SDP to find optimal energy E*.
#   Phase 2: Fix Tr(H ρ) = E* and maximise Tr(ρ²) (purity) to extract
#            the lowest-rank optimizer (extreme point of optimal face).
#            Since Tr(ρ²) is quadratic, we instead add a small random
#            perturbation to break degeneracy and re-solve.
#
# Requires: JuMP, MosekTools, LinearAlgebra, Printf, Random
#

using JuMP, MosekTools, LinearAlgebra, Printf, Random

# ============================================================================
# Hydrogen atom basis for shell n
# ============================================================================

struct HydrogenState
    n::Int
    l::Int
    ml::Int
    ms::Rational{Int}  # ±1//2
end

"""Return all |n, l, mₗ, mₛ⟩ states for the n-th shell (degeneracy 2n²)."""
function shell_basis(n::Int)
    states = HydrogenState[]
    for l in 0:(n-1), ml in -l:l, ms in (-1//2, 1//2)
        push!(states, HydrogenState(n, l, ml, ms))
    end
    return states
end

# ============================================================================
# Operators on a shell
# ============================================================================

function diag_op(basis, f)
    d = length(basis)
    M = zeros(Float64, d, d)
    for (i, s) in enumerate(basis)
        M[i, i] = f(s)
    end
    return M
end

Jz_op(basis) = diag_op(basis, s -> Float64(s.ml + s.ms))
Lz_op(basis) = diag_op(basis, s -> Float64(s.ml))
Sz_op(basis) = diag_op(basis, s -> Float64(s.ms))
L2_op(basis) = diag_op(basis, s -> Float64(s.l * (s.l + 1)))
zeeman_op(basis) = diag_op(basis, s -> Float64(s.ml + 2 * s.ms))

"J² = L² + S² + 2L·S  (includes off-diagonal L₊S₋ + L₋S₊ terms)"
function J2_op(basis)
    d = length(basis)
    M = L2_op(basis) + (3.0/4.0) * I(d)
    M += 2.0 * diag_op(basis, s -> Float64(s.ml) * Float64(s.ms))
    idx = Dict((s.l, s.ml, s.ms) => i for (i, s) in enumerate(basis))
    for (i, s) in enumerate(basis)
        key = (s.l, s.ml + 1, s.ms - 1)
        if haskey(idx, key)
            j = idx[key]
            lp = sqrt(s.l * (s.l + 1) - s.ml * (s.ml + 1))
            sm = sqrt(3/4 - Float64(s.ms) * (Float64(s.ms) - 1))
            M[j, i] += lp * sm
            M[i, j] += lp * sm
        end
    end
    return Hermitian(M)
end

"L·S = (J² - L² - S²)/2"
function fine_structure_op(basis)
    d = length(basis)
    return (Matrix(J2_op(basis)) - Matrix(L2_op(basis)) - (3.0/4.0) * I(d)) / 2.0
end

# ============================================================================
# Two-phase SDP solver
# ============================================================================

"""
    solve_hydrogen_sdp(n, B, constraints; perturbation=1e-8)

Phase 1: min Tr(H ρ) → optimal energy E*.
Phase 2: min Tr((H + εΔ) ρ) with same constraints plus Tr(H ρ) = E*,
         where Δ is a random Hermitian perturbation that generically
         selects a unique extreme point of the optimal face.

Returns (energy, rank, eigenvalues, bp_bound).
"""
function solve_hydrogen_sdp(n::Int, B::Float64,
                            constraints::Vector{<:Tuple};
                            perturbation::Float64=1e-8,
                            seed::Int=42)
    basis = shell_basis(n)
    d = length(basis)
    H = fine_structure_op(basis) + B * zeeman_op(basis)

    # ---- Phase 1: find optimal energy ----
    m1 = Model(Mosek.Optimizer)
    set_silent(m1)
    @variable(m1, ρ1[1:d, 1:d], PSD)
    @constraint(m1, tr(ρ1) == 1.0)
    for (name, A, b) in constraints
        @constraint(m1, tr(A * ρ1) == b)
    end
    @objective(m1, Min, tr(H * ρ1))
    optimize!(m1)

    stat1 = termination_status(m1)
    if stat1 ∉ (OPTIMAL, ALMOST_OPTIMAL, SLOW_PROGRESS)
        return (NaN, -1, Float64[], -1)
    end

    E_star = objective_value(m1)

    # ---- Phase 2: break degeneracy to find extreme point ----
    rng = MersenneTwister(seed)
    Δ = randn(rng, d, d)
    Δ = (Δ + Δ') / 2  # random Hermitian
    H_pert = H + perturbation * Δ

    m2 = Model(Mosek.Optimizer)
    set_silent(m2)
    @variable(m2, ρ2[1:d, 1:d], PSD)
    @constraint(m2, tr(ρ2) == 1.0)
    for (name, A, b) in constraints
        @constraint(m2, tr(A * ρ2) == b)
    end
    # Pin to optimal face
    @constraint(m2, tr(H * ρ2) <= E_star + 1e-7)
    @objective(m2, Min, tr(H_pert * ρ2))
    optimize!(m2)

    stat2 = termination_status(m2)
    if stat2 ∉ (OPTIMAL, ALMOST_OPTIMAL, SLOW_PROGRESS)
        # Fall back to phase 1 result
        ρ_val = Hermitian(value.(ρ1))
    else
        ρ_val = Hermitian(value.(ρ2))
    end

    eigs = sort(real.(eigvals(ρ_val)), rev=true)
    tol = 1e-5
    r = count(e -> e > tol, eigs)
    m_cons = length(constraints)
    bp = floor(Int, sqrt(m_cons + 1))

    return (E_star, r, eigs, bp)
end

# ============================================================================
# Experiment 1: Rank vs number of constraints
# ============================================================================

function experiment_rank_vs_m(n::Int, B::Float64, M::Float64, j_val::Float64, l_val::Int)
    basis = shell_basis(n)
    d = length(basis)
    Jz = Jz_op(basis)
    J2 = Matrix(J2_op(basis))
    L2 = Matrix(L2_op(basis))
    Lz = Lz_op(basis)
    Sz = Sz_op(basis)

    results = NamedTuple[]

    scenarios = [
        (0, Tuple{String,Matrix{Float64},Float64}[], "none"),
        (1, [("Jz", Jz, M)], "Jz"),
        (2, [("Jz", Jz, M), ("J2", J2, j_val*(j_val+1))], "Jz+J2"),
        (3, [("Jz", Jz, M), ("J2", J2, j_val*(j_val+1)),
             ("L2", L2, Float64(l_val*(l_val+1)))], "Jz+J2+L2"),
        (4, [("Jz", Jz, M), ("J2", J2, j_val*(j_val+1)),
             ("L2", L2, Float64(l_val*(l_val+1))),
             ("Lz", Lz, M - 0.5)], "Jz+J2+L2+Lz"),
        (5, [("Jz", Jz, M), ("J2", J2, j_val*(j_val+1)),
             ("L2", L2, Float64(l_val*(l_val+1))),
             ("Lz", Lz, M - 0.5), ("Sz", Sz, 0.5)], "Jz+J2+L2+Lz+Sz"),
    ]

    for (m, cons, label) in scenarios
        E, r, eigs, bp = solve_hydrogen_sdp(n, B, cons)
        if r < 0  # solver failed
            push!(results, (m=m, rank=-1, bp_bound=bp, energy=NaN, scenario=label,
                             eig1=0.0, eig2=0.0, eig3=0.0))
        else
            push!(results, (m=m, rank=r, bp_bound=bp, energy=E, scenario=label,
                             eig1=eigs[1],
                             eig2=length(eigs)>=2 ? eigs[2] : 0.0,
                             eig3=length(eigs)>=3 ? eigs[3] : 0.0))
        end
    end
    return results
end

# ============================================================================
# Experiment 2: Eigenvalue spectrum vs B field
# ============================================================================

function experiment_eigs_vs_B(n::Int, M::Float64, num_constraints::Int,
                               B_range::AbstractVector{Float64},
                               j_val::Float64, l_val::Int)
    basis = shell_basis(n)
    d = length(basis)
    Jz = Jz_op(basis)
    J2 = Matrix(J2_op(basis))
    L2 = Matrix(L2_op(basis))

    results = NamedTuple[]
    for B in B_range
        cons = Tuple{String,Matrix{Float64},Float64}[]
        num_constraints >= 1 && push!(cons, ("Jz", Jz, M))
        num_constraints >= 2 && push!(cons, ("J2", J2, j_val*(j_val+1)))
        num_constraints >= 3 && push!(cons, ("L2", L2, Float64(l_val*(l_val+1))))

        E, r, eigs, bp = solve_hydrogen_sdp(n, B, cons)
        push!(results, (B=B, rank=r, energy=E, bp_bound=bp,
                         eigs=eigs[1:min(d, length(eigs))]))
    end
    return results
end

# ============================================================================
# Experiment 3: Rank vs shell size (dimension scaling)
# ============================================================================

function experiment_rank_vs_shell(B::Float64, M::Float64, m_constraints::Int,
                                   shells::Vector{Int})
    results = NamedTuple[]
    for n in shells
        # Use j=3/2, l=1 when enough constraints; for n=1 shell skip if dim too small
        basis = shell_basis(n)
        d = length(basis)
        if d < 4
            continue
        end

        j_val = 1.5
        l_val = min(1, n-1)

        Jz = Jz_op(basis)
        J2 = Matrix(J2_op(basis))
        L2 = Matrix(L2_op(basis))

        cons = Tuple{String,Matrix{Float64},Float64}[]
        m_constraints >= 1 && push!(cons, ("Jz", Jz, M))
        m_constraints >= 2 && push!(cons, ("J2", J2, j_val*(j_val+1)))
        m_constraints >= 3 && push!(cons, ("L2", L2, Float64(l_val*(l_val+1))))

        E, r, eigs, bp = solve_hydrogen_sdp(n, B, cons)
        push!(results, (shell=n, dim=d, rank=r, bp_bound=bp, energy=E,
                         m=length(cons)))
    end
    return results
end

# ============================================================================
# Data export
# ============================================================================

function export_rank_vs_m(filename, data)
    open(filename, "w") do io
        println(io, "m rank bp_bound energy scenario")
        for r in data
            println(io, @sprintf("%d %d %d %.8f %s", r.m, r.rank, r.bp_bound, r.energy, r.scenario))
        end
    end
end

function export_eigs_vs_B(filename, data, num_eigs=8)
    open(filename, "w") do io
        header = "B rank energy bp_bound " * join(["eig$i" for i in 1:num_eigs], " ")
        println(io, header)
        for r in data
            eigs = vcat(r.eigs, zeros(max(0, num_eigs - length(r.eigs))))
            eig_str = join([@sprintf("%.10f", eigs[i]) for i in 1:num_eigs], " ")
            println(io, @sprintf("%.4f %d %.8f %d %s", r.B, r.rank, r.energy, r.bp_bound, eig_str))
        end
    end
end

function export_rank_vs_shell(filename, data)
    open(filename, "w") do io
        println(io, "shell dim m rank bp_bound energy")
        for r in data
            println(io, @sprintf("%d %d %d %d %d %.8f", r.shell, r.dim, r.m, r.rank, r.bp_bound, r.energy))
        end
    end
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("=" ^ 70)
    println("Barvinok–Pataki rank bounds: Hydrogen atom in magnetic field")
    println("=" ^ 70)

    datadir = joinpath(@__DIR__, "data")
    mkpath(datadir)

    # ---- Experiment 1: rank vs m for M=3/2 (stretched) and M=1/2 ----
    j_val = 1.5
    l_val = 1
    B = 1.0

    println("\n[1/5] Rank vs number of constraints")
    for M in [1.5, 0.5]
        for n in [2, 3, 4]
            println("\n  Shell n=$n (dim=$(2*n^2)), M=$M:")
            println("    m | rank | BP bound | energy      | scenario")
            println("   ---|------|----------|-------------|------------------")
            rv = experiment_rank_vs_m(n, B, M, j_val, l_val)
            for r in rv
                if r.rank < 0
                    println(@sprintf("    %d |  -   |    %d     |         N/A | %s (infeasible)",
                        r.m, r.bp_bound, r.scenario))
                else
                    println(@sprintf("    %d |  %d   |    %d     | %11.6f | %s",
                        r.m, r.rank, r.bp_bound, r.energy, r.scenario))
                end
            end
            tag = @sprintf("M%.0f", 10*M)
            export_rank_vs_m(joinpath(datadir, "rank_vs_m_n$(n)_$(tag).dat"),
                             filter(r -> r.rank >= 0, rv))
        end
    end

    # ---- Experiment 2: eigenvalues of ρ* vs B (n=2 shell) ----
    B_range = collect(0.05:0.05:5.0)
    for M in [1.5, 0.5]
        tag = @sprintf("M%.0f", 10*M)
        println("\n[2/5] Eigenvalue spectrum vs B (n=2, M=$M)")
        for mc in [1, 2, 3]
            evs = experiment_eigs_vs_B(2, M, mc, B_range, j_val, l_val)
            export_eigs_vs_B(joinpath(datadir, "eigs_vs_B_n2_m$(mc)_$(tag).dat"), evs)
            println("  m=$mc: $(length(evs)) data points")
        end
    end

    # ---- Experiment 3: eigenvalues on n=3 shell ----
    println("\n[3/5] Eigenvalue spectrum vs B (n=3, M=0.5)")
    for mc in [1, 3]
        evs = experiment_eigs_vs_B(3, 0.5, mc, B_range, j_val, l_val)
        export_eigs_vs_B(joinpath(datadir, "eigs_vs_B_n3_m$(mc)_M5.dat"), evs, 18)
        println("  m=$mc: $(length(evs)) data points")
    end

    # ---- Experiment 4: rank vs shell size ----
    M = 1.5
    println("\n[4/5] Rank vs shell dimension (M=$M)")
    for mc in [1, 2, 3]
        rv = experiment_rank_vs_shell(B, M, mc, [2, 3, 4, 5])
        export_rank_vs_shell(joinpath(datadir, "rank_vs_shell_m$(mc).dat"), rv)
        println("  m=$mc:")
        for r in rv
            sat = r.rank <= r.bp_bound ? "✓" : "⚠"
            println(@sprintf("    n=%d (dim=%d): rank=%d, BP≤%d %s",
                r.shell, r.dim, r.rank, r.bp_bound, sat))
        end
    end

    # ---- Experiment 5: B sweep for energy and rank (nice for energy plot) ----
    println("\n[5/5] Optimal energy vs B (n=3, m=1,2,3)")
    for mc in [1, 2, 3]
        evs = experiment_eigs_vs_B(3, 1.5, mc, collect(0.1:0.1:5.0), j_val, l_val)
        export_eigs_vs_B(joinpath(datadir, "energy_vs_B_n3_m$(mc).dat"), evs, 18)
    end

    # ---- Summary ----
    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)
    println("Data written to $(datadir)/")
    println("Files:")
    for f in sort(readdir(datadir))
        sz = filesize(joinpath(datadir, f))
        println(@sprintf("  %-35s  %6d bytes", f, sz))
    end
end

main()
