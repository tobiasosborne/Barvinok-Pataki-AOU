#!/usr/bin/env julia
#
# renyi_bp.jl — Rényi entropy at extreme points of spectrahedra
#
# Investigates whether Rényi entropies d_α satisfy bounds tighter than
# the trivial d_α ≤ rank ≤ √m, and how they behave for approximate
# extreme points (robustness).
#
# Three experiments:
#   1. Random spectrahedra: eigenvalue distribution of extreme points
#   2. Robustness: perturb extreme points, track d_α vs perturbation
#   3. Application-specific: channels, marginals — actual d_α vs BP bound
#
# Requires: JuMP, LinearAlgebra, Printf, Random
# Optional: MosekTools (or SCS as fallback)

using JuMP, LinearAlgebra, Printf, Random

# Solver setup
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

# ============================================================================
# Rényi entropy utilities
# ============================================================================

"""Eigenvalues of a density matrix, sorted descending, zeros removed."""
function spectrum(ρ; tol=1e-10)
    λ = eigvals(Hermitian(ρ))
    λ = sort(real.(λ), rev=true)
    return λ[λ .> tol]
end

"""Rényi entropy S_α(ρ) for a density matrix."""
function renyi_entropy(ρ, α; tol=1e-10)
    λ = spectrum(ρ; tol=tol)
    λ = λ / sum(λ)  # renormalise
    if α == 0
        return log(length(λ))
    elseif abs(α - 1) < 1e-12
        return -sum(p * log(p) for p in λ if p > 0)
    elseif α == Inf
        return -log(maximum(λ))
    else
        return log(sum(p^α for p in λ)) / (1 - α)
    end
end

"""Rényi dimension d_α = exp(S_α)."""
renyi_dim(ρ, α; tol=1e-10) = exp(renyi_entropy(ρ, α; tol=tol))

"""Numerical rank."""
numrank(ρ; tol=1e-6) = count(eigvals(Hermitian(ρ)) .> tol)

"""Max r with r² ≤ m (complex BP bound)."""
bp_rank(m) = floor(Int, sqrt(m))

"""Random Hermitian matrix (traceless)."""
function randherm(n)
    M = randn(ComplexF64, n, n)
    H = (M + M') / 2
    return H - tr(H)/n * I
end

"""Random density matrix of given rank."""
function randdm(n, r)
    V = randn(ComplexF64, n, r)
    ρ = V * V'
    return Hermitian(ρ / tr(ρ))
end

# ============================================================================
# Experiment 1: Eigenvalue distribution of extreme points of random spectrahedra
# ============================================================================

"""
Find an extreme point of a random spectrahedron in M_d with m constraints.

Spectrahedron: {ρ ∈ PSD_d : tr(A_i ρ) = b_i, i=1,...,m, tr(ρ)=1}
where A_i are random traceless Hermitian matrices and b_i = tr(A_i ρ₀)
for a random seed state ρ₀ of rank r₀ = bp_rank(m).

Maximises a random linear objective to find an extreme point.
"""
function random_extreme_point(d, m; seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    r₀ = bp_rank(m)
    if r₀ < 1; r₀ = 1; end
    if r₀ > d; r₀ = d; end

    # Generate random constraints
    As = [randherm(d) for _ in 1:m]
    ρ₀ = randdm(d, r₀)
    bs = [real(tr(A * ρ₀)) for A in As]

    # Random objective
    C = randherm(d)

    # Solve SDP
    model = make_model()
    @variable(model, ρ[1:d, 1:d] in HermitianPSDCone())

    for i in 1:m
        @constraint(model, real(tr(As[i] * ρ)) == bs[i])
    end
    @constraint(model, real(tr(ρ)) == 1.0)

    @objective(model, Min, real(tr(C * ρ)))

    optimize!(model)

    if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
        ρ_val = Hermitian(value.(ρ))
        return ρ_val, As, bs
    else
        return nothing, nothing, nothing
    end
end

function experiment1(; d=10, m_values=[4,9,16,25,36,49], n_trials=20)
    println("=" ^ 72)
    println("Experiment 1: Rényi dimensions of extreme points of random spectrahedra")
    println("=" ^ 72)
    println("d = $d, $(n_trials) trials per m value")
    println()

    αs = [0, 0.5, 1, 2, Inf]
    α_names = ["α=0 (rank)", "α=½", "α=1 (vN)", "α=2 (purity)", "α=∞ (1/λ_max)"]

    for m in m_values
        r_bp = bp_rank(m)
        if r_bp > d
            println("m=$m: BP rank $r_bp > d=$d, skipping")
            continue
        end

        println(@sprintf("m = %d,  BP rank bound = %d,  BP d_α bound = %.2f",
                         m, r_bp, sqrt(m)))

        d_α_stats = Dict(α => Float64[] for α in αs)

        for trial in 1:n_trials
            ρ, _, _ = random_extreme_point(d, m; seed=1000*m + trial)
            if ρ === nothing
                continue
            end
            for α in αs
                push!(d_α_stats[α], renyi_dim(ρ, α))
            end
        end

        for (α, name) in zip(αs, α_names)
            vals = d_α_stats[α]
            if isempty(vals); continue; end
            μ = mean(vals)
            σ = length(vals) > 1 ? std(vals) : 0.0
            mx = maximum(vals)
            println(@sprintf("  %-20s  mean=%.3f  std=%.3f  max=%.3f  (bound=%.2f)",
                             name, μ, σ, mx, sqrt(m)))
        end
        println()
    end
end

# ============================================================================
# Experiment 2: Robustness — perturb extreme point, track d_α
# ============================================================================

"""
Perturb an extreme point ρ* by adding noise: ρ_δ = (1-δ)ρ* + δ·I/d.
Track d_α as function of δ.
"""
function experiment2(; d=10, m=16, n_deltas=20)
    println("=" ^ 72)
    println("Experiment 2: Robustness — d_α under depolarising perturbation")
    println("=" ^ 72)

    ρ_star, _, _ = random_extreme_point(d, m; seed=42)
    if ρ_star === nothing
        println("Failed to find extreme point. Skipping.")
        return
    end

    r = numrank(ρ_star)
    println("d=$d, m=$m, rank(ρ*)=$r, BP bound=$(bp_rank(m))")
    println()

    αs = [0, 0.5, 1, 2, Inf]
    δs = [0; 10 .^ range(-8, 0, length=n_deltas)]

    println(@sprintf("%-12s", "δ"), join([@sprintf("%-14s", "d_$(α)") for α in αs]))
    println("-" ^ 82)

    for δ in δs
        ρ_pert = (1 - δ) * ρ_star + δ * I(d) / d
        ρ_pert = Hermitian(ρ_pert / tr(ρ_pert))

        vals = [renyi_dim(ρ_pert, α; tol=1e-14) for α in αs]
        println(@sprintf("%-12.2e", δ),
                join([@sprintf("%-14.4f", v) for v in vals]))
    end

    println()
    println("Note: d_0 (rank) jumps from $r to $d at any δ>0.")
    println("      d_α for α>0 changes continuously.")
    println("      BP bound: d_α ≤ $(bp_rank(m)) for all α at the exact extreme point.")
end

# ============================================================================
# Experiment 3: Channels — Rényi entropies of extreme CPTP maps
# ============================================================================

"""
Find an extreme CPTP map E: M_{d_in} → M_{d_out} by optimising over
Choi matrices with CPTP constraints.

CPTP constraints on Choi matrix J ∈ PSD_{d_out·d_in}:
  J ≥ 0,  tr_{d_out}(J) = I_{d_in}
The number of real constraints is m = d_in².
BP bound on rank(J): r ≤ d_in (= √m).
"""
function extreme_channel(d_in, d_out; seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    d = d_out * d_in

    # Random objective
    C = randherm(d)

    model = make_model()
    @variable(model, J[1:d, 1:d] in HermitianPSDCone())

    # Partial trace over output = I_in
    for a in 1:d_in, b in 1:d_in
        val = (a == b) ? 1.0 : 0.0
        # tr_{d_out}(J)_{a,b} = Σ_k J[(k-1)*d_in+a, (k-1)*d_in+b]
        expr = sum(J[(k-1)*d_in+a, (k-1)*d_in+b] for k in 1:d_out)
        @constraint(model, real(expr) == val)
        if a != b
            @constraint(model, imag(expr) == 0.0)
        end
    end

    @objective(model, Min, real(tr(C * J)))
    optimize!(model)

    if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
        J_val = Hermitian(value.(J))
        # Normalise Choi matrix to be a state: ρ_J = J / d_in
        ρ_J = J_val / d_in
        return ρ_J, J_val
    else
        return nothing, nothing
    end
end

function experiment3(; dims=[(2,2), (2,3), (2,4), (3,3), (3,4)], n_trials=10)
    println("=" ^ 72)
    println("Experiment 3: Rényi entropies of extreme CPTP channels")
    println("=" ^ 72)
    println()

    αs = [0, 1, 2, Inf]

    for (d_in, d_out) in dims
        m = d_in^2  # number of CPTP constraints
        r_bp = d_in  # BP bound on Kraus rank
        d_choi = d_in * d_out

        println(@sprintf("d_in=%d, d_out=%d  |  Choi dim=%d, m=%d, BP Kraus rank ≤ %d",
                         d_in, d_out, d_choi, m, r_bp))

        for trial in 1:n_trials
            ρ_J, _ = extreme_channel(d_in, d_out; seed=100*d_in+10*d_out+trial)
            if ρ_J === nothing; continue; end

            r = numrank(ρ_J)
            d_vals = [renyi_dim(ρ_J, α) for α in αs]

            if trial == 1
                println(@sprintf("  Trial %2d: rank=%d  d₀=%.1f  d₁=%.3f  d₂=%.3f  d_∞=%.3f",
                                 trial, r, d_vals...))
            elseif trial == n_trials
                println(@sprintf("  Trial %2d: rank=%d  d₀=%.1f  d₁=%.3f  d₂=%.3f  d_∞=%.3f",
                                 trial, r, d_vals...))
            end
        end
        println()
    end
end

# ============================================================================
# Experiment 4: Eigenvalue distribution at saturation
# ============================================================================

"""
For extreme points that saturate BP (rank = √m), are eigenvalues
typically uniform or skewed? This determines whether d_α < d₀.
"""
function experiment4(; d=12, m=16, n_trials=50)
    println("=" ^ 72)
    println("Experiment 4: Eigenvalue skewness at BP-saturating extreme points")
    println("=" ^ 72)

    r_bp = bp_rank(m)
    println("d=$d, m=$m, BP rank bound=$r_bp")
    println()

    n_saturating = 0
    skewness_data = Float64[]  # ratio d_2/d_0 (1 = uniform, <1 = skewed)

    for trial in 1:n_trials
        ρ, _, _ = random_extreme_point(d, m; seed=2000+trial)
        if ρ === nothing; continue; end

        r = numrank(ρ)
        if r == r_bp
            n_saturating += 1
            ratio = renyi_dim(ρ, 2) / renyi_dim(ρ, 0)
            push!(skewness_data, ratio)

            λ = spectrum(ρ)
            if n_saturating <= 5
                println(@sprintf("  Saturating trial: rank=%d, eigenvalues=[%s]",
                                 r, join([@sprintf("%.4f", l) for l in λ], ", ")))
                println(@sprintf("    d₂/d₀ = %.4f  (1.0 = uniform)", ratio))
            end
        end
    end

    println()
    println(@sprintf("BP-saturating: %d / %d trials", n_saturating, n_trials))
    if !isempty(skewness_data)
        println(@sprintf("d₂/d₀ ratio:  mean=%.4f  std=%.4f  min=%.4f  max=%.4f",
                         mean(skewness_data), std(skewness_data),
                         minimum(skewness_data), maximum(skewness_data)))
        println()
        if mean(skewness_data) < 0.95
            println(">>> Eigenvalues are typically SKEWED (not uniform).")
            println("    This means d_α < d₀ = √m generically,")
            println("    so the Rényi bound IS tighter than the rank bound in practice!")
        else
            println(">>> Eigenvalues are typically UNIFORM.")
            println("    The Rényi bound equals the rank bound generically.")
        end
    end
end

# ============================================================================
# Experiment 5: Pataki eigenvalue clustering at SDP optima
# ============================================================================

"""
Test whether optimising a linear objective min tr(C·ρ) over a spectrahedron
produces eigenvalue clustering at the optimum compared to a random feasible
extreme point.  Measures eigenvalue gaps and d₂/d₀ ratios.
"""
function experiment5(; d=12, m=16, n_trials=30, gap_tol=1e-4)
    println("=" ^ 72)
    println("Experiment 5: Pataki eigenvalue clustering at SDP optima")
    println("=" ^ 72)
    println("d=$d, m=$m, gap tolerance=$gap_tol, $n_trials trials")
    println()

    opt_clusters   = Int[]
    rand_clusters  = Int[]
    opt_ratios     = Float64[]
    rand_ratios    = Float64[]

    for trial in 1:n_trials
        # Optimised extreme point (primary objective)
        ρ_opt, As, bs = random_extreme_point(d, m; seed=5000 + trial)
        if ρ_opt === nothing; continue; end

        # Random extreme point (different objective, same spectrahedron)
        Random.seed!(6000 + trial)
        C2 = randherm(d)
        model2 = make_model()
        @variable(model2, ρ2[1:d, 1:d] in HermitianPSDCone())
        for i in 1:m
            @constraint(model2, real(tr(As[i] * ρ2)) == bs[i])
        end
        @constraint(model2, real(tr(ρ2)) == 1.0)
        @objective(model2, Min, real(tr(C2 * ρ2)))
        optimize!(model2)

        if !(termination_status(model2) in [OPTIMAL, ALMOST_OPTIMAL]); continue; end
        ρ_rand = Hermitian(value.(ρ2))

        # Eigenvalue gaps
        λ_opt  = spectrum(ρ_opt)
        λ_rand = spectrum(ρ_rand)

        gaps_opt  = [abs(λ_opt[i] - λ_opt[i+1]) for i in 1:length(λ_opt)-1]
        gaps_rand = [abs(λ_rand[i] - λ_rand[i+1]) for i in 1:length(λ_rand)-1]

        n_cluster_opt  = count(g -> g < gap_tol, gaps_opt)
        n_cluster_rand = count(g -> g < gap_tol, gaps_rand)

        push!(opt_clusters,  n_cluster_opt)
        push!(rand_clusters, n_cluster_rand)

        # d₂/d₀ ratio
        r_opt  = renyi_dim(ρ_opt, 2) / renyi_dim(ρ_opt, 0)
        r_rand = renyi_dim(ρ_rand, 2) / renyi_dim(ρ_rand, 0)
        push!(opt_ratios,  r_opt)
        push!(rand_ratios, r_rand)

        if trial <= 5
            println(@sprintf("  Trial %2d: opt  rank=%d  clusters=%d  d₂/d₀=%.4f",
                             trial, numrank(ρ_opt), n_cluster_opt, r_opt))
            println(@sprintf("            rand rank=%d  clusters=%d  d₂/d₀=%.4f",
                             numrank(ρ_rand), n_cluster_rand, r_rand))
        end
    end

    println()
    n = length(opt_clusters)
    if n > 0
        println(@sprintf("Summary over %d successful trials:", n))
        println(@sprintf("  Optimum:  mean clusters=%.2f  mean d₂/d₀=%.4f ± %.4f",
                         mean(opt_clusters), mean(opt_ratios), std(opt_ratios)))
        println(@sprintf("  Random:   mean clusters=%.2f  mean d₂/d₀=%.4f ± %.4f",
                         mean(rand_clusters), mean(rand_ratios), std(rand_ratios)))
        println()
        if mean(opt_clusters) > mean(rand_clusters) + 0.5
            println(">>> Optima show MORE eigenvalue clustering than random extreme points.")
        else
            println(">>> No significant difference in clustering between optima and random.")
        end
    end
end

# ============================================================================
# Experiment 6: Scaling of d_α with m (test typical bound conjecture)
# ============================================================================

"""
For fixed d and varying m, find extreme points and compute d_α / √m for
α ∈ {1, 2, ∞}.  Tests whether d_α ≤ c_α · √m with c_α < 1.
"""
function experiment6(; d=20, m_values=[4, 9, 16, 25, 36, 49, 64, 81, 100], n_trials=10)
    println("=" ^ 72)
    println("Experiment 6: Scaling of d_α with m — typical bound conjecture")
    println("=" ^ 72)
    println("d=$d, $n_trials trials per m")
    println()

    αs = [1, 2, Inf]
    α_labels = ["c₁ (vN)", "c₂ (purity)", "c_∞"]

    println(@sprintf("%-8s  %-14s  %-14s  %-14s  %-10s",
                     "m", α_labels[1], α_labels[2], α_labels[3], "√m"))
    println("-" ^ 66)

    for m in m_values
        r_bp = bp_rank(m)
        if r_bp > d
            println(@sprintf("%-8d  (BP rank %d > d=%d, skipped)", m, r_bp, d))
            continue
        end

        ratios = Dict(α => Float64[] for α in αs)

        for trial in 1:n_trials
            ρ, _, _ = random_extreme_point(d, m; seed=7000 + 100*m + trial)
            if ρ === nothing; continue; end

            for α in αs
                dα = renyi_dim(ρ, α)
                push!(ratios[α], dα / sqrt(m))
            end
        end

        c_vals = [isempty(ratios[α]) ? NaN : mean(ratios[α]) for α in αs]
        println(@sprintf("%-8d  %-14.4f  %-14.4f  %-14.4f  %-10.2f",
                         m, c_vals..., sqrt(m)))
    end

    println()
    println("If c_α < 1 consistently, then d_α < √m generically (tighter than BP rank).")
end

# ============================================================================
# Experiment 7: Application-specific Rényi dimensions
# ============================================================================

"""
Compute d_α for extreme points of:
  a) Marginal-constrained bipartite states
  b) Covariant qubit channels
"""
function experiment7(; d_A=3, d_B=3)
    println("=" ^ 72)
    println("Experiment 7: Application-specific Rényi dimensions")
    println("=" ^ 72)
    println()

    # --- Part (a): Marginal-constrained bipartite states ---
    println("--- (a) Marginal-constrained bipartite states ---")
    println("d_A=$d_A, d_B=$d_B")

    d = d_A * d_B
    m = d_A^2 + d_B^2 - 1  # constraints from both marginals (trace auto-satisfied)

    Random.seed!(8001)
    σ_A = randdm(d_A, d_A)
    σ_B = randdm(d_B, d_B)

    C = randherm(d)

    model = make_model()
    @variable(model, ρ[1:d, 1:d] in HermitianPSDCone())

    # Constraint: tr_B(ρ) = σ_A
    for a in 1:d_A, b in 1:d_A
        # tr_B(ρ)_{a,b} = Σ_k ρ[(a-1)*d_B+k, (b-1)*d_B+k]
        expr = sum(ρ[(a-1)*d_B+k, (b-1)*d_B+k] for k in 1:d_B)
        @constraint(model, real(expr) == real(σ_A[a,b]))
        if a != b
            @constraint(model, imag(expr) == imag(σ_A[a,b]))
        end
    end

    # Constraint: tr_A(ρ) = σ_B
    for a in 1:d_B, b in 1:d_B
        # tr_A(ρ)_{a,b} = Σ_k ρ[(k-1)*d_B+a, (k-1)*d_B+b]
        expr = sum(ρ[(k-1)*d_B+a, (k-1)*d_B+b] for k in 1:d_A)
        @constraint(model, real(expr) == real(σ_B[a,b]))
        if a != b
            @constraint(model, imag(expr) == imag(σ_B[a,b]))
        end
    end

    @objective(model, Min, real(tr(C * ρ)))
    optimize!(model)

    if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
        ρ_val = Hermitian(value.(ρ))
        r = numrank(ρ_val)
        d1 = renyi_dim(ρ_val, 1)
        d2 = renyi_dim(ρ_val, 2)
        d_inf = renyi_dim(ρ_val, Inf)
        r_bp = bp_rank(m)
        println(@sprintf("  m=%d constraints, BP rank ≤ %d", m, r_bp))
        println(@sprintf("  rank=%d,  d₁=%.4f,  d₂=%.4f,  d_∞=%.4f", r, d1, d2, d_inf))
        println(@sprintf("  BP bound: d_α ≤ √m = %.4f", sqrt(m)))
    else
        println("  Marginal SDP failed: $(termination_status(model))")
    end
    println()

    # --- Part (b): Covariant qubit channels ---
    println("--- (b) Covariant qubit channels (σ_z symmetry) ---")

    d_in = 2
    d_out = 2
    d_choi = d_in * d_out  # = 4

    Random.seed!(8002)
    C_ch = randherm(d_choi)

    model2 = make_model()
    @variable(model2, J[1:d_choi, 1:d_choi] in HermitianPSDCone())

    # CPTP constraints: tr_{out}(J) = I_{in}
    for a in 1:d_in, b in 1:d_in
        val = (a == b) ? 1.0 : 0.0
        expr = sum(J[(k-1)*d_in+a, (k-1)*d_in+b] for k in 1:d_out)
        @constraint(model2, real(expr) == val)
        if a != b
            @constraint(model2, imag(expr) == 0.0)
        end
    end

    # Covariance constraint: [E(·), σ_z] = 0
    # On the Choi matrix: (σ_z ⊗ I) J = J (σ_z ⊗ I)
    # σ_z = diag(1, -1)
    σ_z = [1.0 0.0; 0.0 -1.0]
    σ_z_I = kron(σ_z, Matrix{Float64}(I, d_in, d_in))

    comm = σ_z_I  # (σ_z ⊗ I)
    for i in 1:d_choi, j in 1:d_choi
        # [comm, J]_{i,j} = Σ_k comm_{i,k} J_{k,j} - J_{i,k} comm_{k,j}
        expr = sum(comm[i,k] * J[k,j] - J[i,k] * comm[k,j] for k in 1:d_choi)
        @constraint(model2, real(expr) == 0.0)
        @constraint(model2, imag(expr) == 0.0)
    end

    @objective(model2, Min, real(tr(C_ch * J)))
    optimize!(model2)

    if termination_status(model2) in [OPTIMAL, ALMOST_OPTIMAL]
        J_val = Hermitian(value.(J))
        ρ_J = J_val / d_in
        r = numrank(ρ_J)
        d1 = renyi_dim(ρ_J, 1)
        d2 = renyi_dim(ρ_J, 2)
        d_inf = renyi_dim(ρ_J, Inf)
        m_cptp = d_in^2
        m_cov = 2 * d_choi^2  # real + imag parts of commutator
        m_total = m_cptp + m_cov
        r_bp = bp_rank(m_total)
        println(@sprintf("  d_in=%d, d_out=%d, Choi dim=%d", d_in, d_out, d_choi))
        println(@sprintf("  CPTP constraints=%d, covariance constraints=%d, total m=%d",
                         m_cptp, m_cov, m_total))
        println(@sprintf("  rank=%d,  d₁=%.4f,  d₂=%.4f,  d_∞=%.4f", r, d1, d2, d_inf))
        println(@sprintf("  BP bound: d_α ≤ √m = %.4f", sqrt(m_total)))
    else
        println("  Covariant channel SDP failed: $(termination_status(model2))")
    end
end

# ============================================================================
# Helpers
# ============================================================================
mean(x) = sum(x) / length(x)
function std(x)
    μ = mean(x)
    return sqrt(sum((xi - μ)^2 for xi in x) / max(length(x) - 1, 1))
end

# ============================================================================
# Run all experiments
# ============================================================================

function main()
    println("Rényi Generalisations of the Barvinok–Pataki Bound")
    println("Numerical Investigation — $(Dates.now())")
    println()

    experiment1(d=10, m_values=[4, 9, 16, 25], n_trials=15)
    println()
    experiment2(d=10, m=16)
    println()
    experiment3(dims=[(2,2), (2,3), (3,3)], n_trials=5)
    println()
    experiment4(d=12, m=16, n_trials=30)
    println()
    experiment5(d=12, m=16, n_trials=30)
    println()
    experiment6(d=20, m_values=[4, 9, 16, 25, 36, 49, 64, 81, 100], n_trials=10)
    println()
    experiment7(d_A=3, d_B=3)
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Dates
    main()
end
