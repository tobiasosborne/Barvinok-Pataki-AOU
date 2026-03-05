#!/usr/bin/env julia
#
# scattering_bp.jl — Barvinok–Pataki for inclusive scattering channels
#
# Three scattering examples:
#   1. Generic 2→2 with spin: parametric in d_in, k conservation laws
#   2. Electron-electron (Coulomb): partial-wave truncation, J_z + parity
#   3. Compton scattering: γ+e⁻→γ+e⁻, helicity basis, trace out photon
#
# The S-matrix is the Stinespring dilation of the inclusive channel
# Φ_S(ρ) = Tr_unobs(S ρ S†). Conservation laws are linear constraints
# on the Choi matrix, tightening the BP bound on Kraus rank.
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

"""Max r such that r² ≤ m (complex Hermitian BP bound)."""
bp_rank(m) = floor(Int, sqrt(m))

"""Max r such that r(r+1)/2 ≤ m (real symmetric BP bound)."""
bp_rank_real(m) = floor(Int, (-1 + sqrt(1 + 8m)) / 2)

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


# ============================================================================
# Example 1: Generic 2→2 scattering with spin
# ============================================================================
#
# Two incoming particles → two outgoing. Detect one outgoing particle.
# d_in = input Hilbert space dimension (momentum × spin)
# d_det = detected particle dimension
# d_unobs = unobserved particle dimension (d_out = d_det · d_unobs)
#
# Inclusive channel: Φ(ρ) = Tr_unobs(S ρ S†)
# Choi matrix J ∈ PSD_{d_in · d_det} with Tr_{d_det}(J) = I_{d_in}
# → m = d_in² constraints (TP condition)
#
# k conservation laws: Tr(Q_i J) = q_i → m_total = d_in² + k
# BP: r ≤ √(d_in² + k)
#
function example1_generic_scattering()
    println("\n" * "="^72)
    println("Example 1: Generic 2→2 scattering — parametric BP bounds")
    println("="^72)

    println("\n  Table: Kraus rank bound for inclusive channel")
    println("  " * "-"^60)
    println("  d_in  d_det | k=0  k=1  k=2  k=3  k=5  k=10")
    println("  " * "-"^60)

    for d_in in [2, 3, 4, 6, 8, 10]
        d_det = d_in  # same dimension for detected particle
        line = @sprintf("  %3d   %3d   |", d_in, d_det)
        for k in [0, 1, 2, 3, 5, 10]
            m_total = d_in^2 + k
            r_max = bp_rank(m_total)
            line *= @sprintf(" %3d ", r_max)
        end
        println(line)
    end
    println("  " * "-"^60)
    println("  (r ≤ ⌊√(d_in² + k)⌋)")

    # Numerical verification for d_in = 3
    println("\n  Numerical verification (d_in=3, d_det=3):")
    Random.seed!(42)
    dA = 3; dB = 3; d = dA * dB

    for k in [0, 1, 2, 3]
        m_tp = dA^2
        m_total = m_tp + k
        bp_bound = bp_rank(m_total)

        max_rank = 0
        for trial in 1:8
            model = make_model()
            @variable(model, J[1:d, 1:d], PSD)

            # TP: Tr_B(J) = I_A
            for ia in 1:dA, ja in ia:dA
                lhs = @expression(model, sum(J[(ia-1)*dB+ib, (ja-1)*dB+ib] for ib in 1:dB))
                rhs = (ia == ja) ? 1.0 : 0.0
                @constraint(model, lhs == rhs)
            end

            # k conservation law constraints (random Hermitian operators on d_in ⊗ d_det)
            for c in 1:k
                Q = let M = randn(d, d); (M + M') / 2 end
                q = tr(Q) / d  # value at maximally mixed
                @constraint(model, sum(Q[i,j] * J[i,j] for i in 1:d, j in 1:d) == q)
            end

            C = let M = randn(d, d); (M + M') / 2 end
            @objective(model, Min, sum(C[i,j] * J[i,j] for i in 1:d, j in 1:d))

            optimize!(model)
            if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
                Jval = value.(J); Jval = (Jval + Jval') / 2
                max_rank = max(max_rank, numrank(Jval))
            end
        end
        @printf("    k=%d: m=%d, BP r≤%d, observed max rank=%d  %s\n",
                k, m_total, bp_bound, max_rank, max_rank <= bp_bound ? "✓" : "✗")
    end
end


# ============================================================================
# Example 2: Electron-electron Coulomb scattering
# ============================================================================
#
# e⁻ + e⁻ → e⁻ + e⁻ (identical fermions)
# Partial-wave decomposition: truncate at L ≤ L_max
# Spin: each electron is spin-1/2
#   Total spin S ∈ {0, 1} (singlet/triplet)
# Angular momentum: L partial waves coupled with spin
#   |L, M_L⟩ ⊗ |S, M_S⟩, but we work with J, J_z basis:
#   |J, J_z⟩ with J = |L-S|,...,L+S
#
# Conservation laws:
#   - Total angular momentum J_z (always conserved)
#   - Parity P = (-1)^L (always conserved for Coulomb)
#   - Optionally: total J (if we enforce it)
#
# Inclusive channel: trace out one electron (detected vs unobserved)
#
function example2_electron_electron()
    println("\n" * "="^72)
    println("Example 2: Electron-electron Coulomb scattering")
    println("="^72)

    println("\n  Partial-wave truncation, J_z and parity conservation")
    println("  " * "-"^60)
    println("  L_max | d_in | d_det | m_TP  | m_Jz | m_P | m_tot | r_BP | rank_obs")
    println("  " * "-"^60)

    Random.seed!(42)

    for Lmax in [1, 2, 3, 4]
        # Build Hilbert space in |L, M_L, S, M_S⟩ basis
        # For identical fermions: L+S must be even (antisymmetry)
        # Relative motion states: |L, M_L⟩ with L = 0,...,L_max
        # Spin states: |S, M_S⟩ with S=0 (singlet) or S=1 (triplet)

        basis_labels = []  # (L, M_L, S, M_S)
        for L in 0:Lmax
            for ML in -L:L
                for S in [0, 1]
                    if (L + S) % 2 != 0; continue; end  # antisymmetry
                    for MS in -S:S
                        push!(basis_labels, (L, ML, S, MS))
                    end
                end
            end
        end

        d_full = length(basis_labels)
        if d_full == 0; continue; end

        # For 2→2 scattering of identical particles, input = output space
        # Inclusive: detect one electron → need to split into "detected" and "unobserved"
        # In partial-wave basis, the S-matrix is block-diagonal in (J, J_z, P).
        # The inclusive channel traces out the relative-motion quantum numbers
        # while keeping the detected electron's spin.
        #
        # Simplification: treat the full partial-wave space as the input,
        # and the detected particle as a 2-dimensional spin-1/2 system.
        dA = d_full   # input (full scattering space)
        dB = 2        # detected electron spin
        d = dA * dB

        m_tp = dA * (dA + 1) ÷ 2   # TP constraints (real symmetric Choi)
        # In practice we cap at reasonable sizes
        if d > 40
            # Compute constraints analytically — iterate over Choi upper triangle
            # exactly as the SDP branch does, but without building the model
            Jz_vals = [basis_labels[i][2] + basis_labels[i][4] for i in 1:d_full]  # ML + MS
            P_vals = [(-1.0)^basis_labels[i][1] for i in 1:d_full]  # (-1)^L
            Sz_vals = [0.5, -0.5]  # detected spin-z

            n_Jz = 0
            for ia in 1:dA, a in 1:dB, ja in 1:dA, b in 1:dB
                i_idx = (ia-1)*dB + a; j_idx = (ja-1)*dB + b
                if i_idx > j_idx; continue; end
                Jz_i = Jz_vals[ia] + Sz_vals[a]
                Jz_j = Jz_vals[ja] + Sz_vals[b]
                if abs(Jz_i - Jz_j) > 1e-10
                    n_Jz += 1
                end
            end

            n_P = 0
            for ia in 1:dA, a in 1:dB, ja in 1:dA, b in 1:dB
                i_idx = (ia-1)*dB + a; j_idx = (ja-1)*dB + b
                if i_idx > j_idx; continue; end
                if abs(P_vals[ia] - P_vals[ja]) > 1e-10
                    # Skip if already constrained by J_z
                    Jz_i = Jz_vals[ia] + Sz_vals[a]
                    Jz_j = Jz_vals[ja] + Sz_vals[b]
                    if abs(Jz_i - Jz_j) > 1e-10; continue; end
                    n_P += 1
                end
            end

            m_total = m_tp + n_Jz + n_P
            r_bp = bp_rank_real(m_total)
            @printf("  %3d   | %3d  | %3d   | %5d | %4d | %3d | %5d | %4d | (skip: d=%d)\n",
                    Lmax, dA, dB, m_tp, n_Jz, n_P, m_total, r_bp, d)
            continue
        end

        # Build conservation law operators on the Choi space
        # J_z operator on input space
        Jz_in = zeros(dA, dA)
        for i in 1:dA
            L, ML, S, MS = basis_labels[i]
            Jz_in[i, i] = ML + MS
        end

        # S_z operator on detected spin-1/2
        Sz_det = [0.5 0.0; 0.0 -0.5]

        # J_z on Choi space: J_z^in ⊗ I_det + I_in ⊗ S_z^det
        # Conservation: [J_z_total, J_Choi] = 0
        # This means J_Choi[(i,a),(j,b)] = 0 if (Jz_in[i]+Sz[a]) ≠ (Jz_in[j]+Sz[b])

        # Parity operator on input space
        P_in = zeros(dA, dA)
        for i in 1:dA
            L = basis_labels[i][1]
            P_in[i, i] = (-1.0)^L
        end

        # Build SDP
        model = make_model()
        @variable(model, J_choi[1:d, 1:d], PSD)

        # TP: Tr_B(J) = I_A
        for ia in 1:dA, ja in ia:dA
            lhs = @expression(model, sum(J_choi[(ia-1)*dB+ib, (ja-1)*dB+ib] for ib in 1:dB))
            rhs = (ia == ja) ? 1.0 : 0.0
            @constraint(model, lhs == rhs)
        end

        # J_z conservation constraints
        n_Jz = 0
        Sz_vals = [0.5, -0.5]
        for ia in 1:dA, a in 1:dB, ja in 1:dA, b in 1:dB
            i = (ia-1)*dB + a
            j = (ja-1)*dB + b
            if i > j; continue; end
            Jz_i = Jz_in[ia, ia] + Sz_vals[a]
            Jz_j = Jz_in[ja, ja] + Sz_vals[b]
            if abs(Jz_i - Jz_j) > 1e-10
                @constraint(model, J_choi[i, j] == 0.0)
                n_Jz += 1
            end
        end

        # Parity conservation constraints
        n_P = 0
        for ia in 1:dA, a in 1:dB, ja in 1:dA, b in 1:dB
            i = (ia-1)*dB + a
            j = (ja-1)*dB + b
            if i > j; continue; end
            Pi = P_in[ia, ia]
            Pj = P_in[ja, ja]
            if abs(Pi - Pj) > 1e-10
                # Already constrained by J_z? Check and skip if so
                Jz_i = Jz_in[ia, ia] + Sz_vals[a]
                Jz_j = Jz_in[ja, ja] + Sz_vals[b]
                if abs(Jz_i - Jz_j) > 1e-10
                    continue  # already zero from J_z
                end
                @constraint(model, J_choi[i, j] == 0.0)
                n_P += 1
            end
        end

        m_total = dA * (dA + 1) ÷ 2 + n_Jz + n_P  # upper-triangle TP + conservation
        bp_bound = bp_rank_real(m_total)

        # Solve with random objectives
        max_rank = 0
        for trial in 1:8
            model_t = make_model()
            @variable(model_t, Jt[1:d, 1:d], PSD)

            for ia in 1:dA, ja in ia:dA
                lhs = @expression(model_t, sum(Jt[(ia-1)*dB+ib, (ja-1)*dB+ib] for ib in 1:dB))
                rhs = (ia == ja) ? 1.0 : 0.0
                @constraint(model_t, lhs == rhs)
            end
            for ia in 1:dA, a in 1:dB, ja in 1:dA, b in 1:dB
                i = (ia-1)*dB + a; j = (ja-1)*dB + b
                if i > j; continue; end
                Jz_i = Jz_in[ia, ia] + Sz_vals[a]
                Jz_j = Jz_in[ja, ja] + Sz_vals[b]
                if abs(Jz_i - Jz_j) > 1e-10
                    @constraint(model_t, Jt[i, j] == 0.0)
                end
            end
            for ia in 1:dA, a in 1:dB, ja in 1:dA, b in 1:dB
                i = (ia-1)*dB + a; j = (ja-1)*dB + b
                if i > j; continue; end
                Pi = P_in[ia, ia]; Pj = P_in[ja, ja]
                if abs(Pi - Pj) > 1e-10
                    Jz_i = Jz_in[ia, ia] + Sz_vals[a]
                    Jz_j = Jz_in[ja, ja] + Sz_vals[b]
                    if abs(Jz_i - Jz_j) > 1e-10; continue; end
                    @constraint(model_t, Jt[i, j] == 0.0)
                end
            end

            C = let M = randn(d, d); (M + M') / 2 end
            @objective(model_t, Min, sum(C[i,j] * Jt[i,j] for i in 1:d, j in 1:d))
            optimize!(model_t)
            if termination_status(model_t) in [OPTIMAL, ALMOST_OPTIMAL]
                Jval = value.(Jt); Jval = (Jval + Jval') / 2
                max_rank = max(max_rank, numrank(Jval))
            end
        end

        @printf("  %3d   | %3d  | %3d   | %5d | %4d | %3d | %5d | %4d | %4d  %s\n",
                Lmax, dA, dB, dA*(dA+1)÷2, n_Jz, n_P, m_total, bp_bound,
                max_rank, max_rank <= bp_bound ? "✓" : "✗")
    end
end


# ============================================================================
# Example 3: Compton scattering (γ + e⁻ → γ + e⁻)
# ============================================================================
#
# Photon: spin-1 (helicity ±1, so effective d=2)
# Electron: spin-1/2 (d=2)
# Partial-wave decomposition in helicity basis
#
# Inclusive: detect electron, trace out photon
# Channel Φ: M_{d_in} → M_{d_electron} (d_electron=2)
#
# Conservation: helicity (J_z), parity
#
function example3_compton()
    println("\n" * "="^72)
    println("Example 3: Compton scattering (γ + e⁻ → γ + e⁻)")
    println("="^72)

    println("\n  Helicity basis, partial-wave truncation, detect electron")
    println("  " * "-"^60)
    println("  J_max | d_in | d_det | m_TP  | m_cons | m_tot | r_BP | rank_obs")
    println("  " * "-"^60)

    Random.seed!(77)

    for Jmax in [1, 2, 3, 4]
        # Helicity basis for γ + e⁻:
        # Photon helicities: λ_γ ∈ {+1, -1}
        # Electron helicities: λ_e ∈ {+1/2, -1/2}
        # Total helicity: λ = λ_γ + λ_e
        # Partial waves: J ≥ |λ|, J = 1/2, 3/2, ..., Jmax (half-integer)
        # Basis: |J, λ⟩ for each valid (J, λ_γ, λ_e) combination

        # For the scattering problem, we label states by |J, λ_γ, λ_e⟩
        # where J ≥ |λ_γ + λ_e| and J is half-integer
        basis = []  # (J, λ_γ, λ_e)
        for λγ in [1, -1]
            for λe in [0.5, -0.5]
                λ_total = λγ + λe
                # J must be half-integer ≥ |λ_total|
                J_min = abs(λ_total)
                if J_min < 0.5; J_min = 0.5; end
                # Make sure J_min is half-integer
                J_min = ceil(J_min - 0.5) + 0.5
                J = J_min
                while J <= Jmax + 0.5
                    push!(basis, (J, λγ, λe))
                    J += 1.0
                end
            end
        end

        d_in = length(basis)
        if d_in == 0; continue; end
        d_det = 2  # detected electron spin
        d = d_in * d_det

        if d > 50
            # Analytical bounds only — iterate over Choi upper triangle
            m_tp = d_in * (d_in + 1) ÷ 2
            λ_in = [basis[i][2] + basis[i][3] for i in 1:d_in]
            Sz_vals = [0.5, -0.5]
            n_cons = 0
            for ia in 1:d_in, a in 1:d_det, ja in 1:d_in, b in 1:d_det
                i_idx = (ia-1)*d_det + a; j_idx = (ja-1)*d_det + b
                if i_idx > j_idx; continue; end
                λ_i = λ_in[ia] + Sz_vals[a]
                λ_j = λ_in[ja] + Sz_vals[b]
                if abs(λ_i - λ_j) > 1e-10
                    n_cons += 1
                end
            end
            m_total = m_tp + n_cons
            r_bp = bp_rank_real(m_total)
            @printf("  %3.0f   | %3d  | %3d   | %5d | %5d  | %5d | %4d | (skip: d=%d)\n",
                    Jmax, d_in, d_det, m_tp, n_cons, m_total, r_bp, d)
            continue
        end

        # Build J_z (helicity) on input space
        Jz_in = zeros(d_in, d_in)
        for i in 1:d_in
            J_val, λγ, λe = basis[i]
            Jz_in[i, i] = λγ + λe  # total helicity
        end

        Sz_vals = [0.5, -0.5]  # detected electron spin-z

        # SDP
        max_rank = 0
        n_Jz = 0  # count conservation constraints

        for trial in 1:8
            model = make_model()
            @variable(model, Jc[1:d, 1:d], PSD)

            # TP
            for ia in 1:d_in, ja in ia:d_in
                lhs = @expression(model, sum(Jc[(ia-1)*d_det+ib, (ja-1)*d_det+ib] for ib in 1:d_det))
                rhs = (ia == ja) ? 1.0 : 0.0
                @constraint(model, lhs == rhs)
            end

            # Helicity conservation
            nc = 0
            for ia in 1:d_in, a in 1:d_det, ja in 1:d_in, b in 1:d_det
                i = (ia-1)*d_det + a; j = (ja-1)*d_det + b
                if i > j; continue; end
                λ_i = Jz_in[ia, ia] + Sz_vals[a]
                λ_j = Jz_in[ja, ja] + Sz_vals[b]
                if abs(λ_i - λ_j) > 1e-10
                    @constraint(model, Jc[i, j] == 0.0)
                    nc += 1
                end
            end
            if trial == 1; n_Jz = nc; end

            C = let M = randn(d, d); (M + M') / 2 end
            @objective(model, Min, sum(C[i,j] * Jc[i,j] for i in 1:d, j in 1:d))
            optimize!(model)
            if termination_status(model) in [OPTIMAL, ALMOST_OPTIMAL]
                Jval = value.(Jc); Jval = (Jval + Jval') / 2
                max_rank = max(max_rank, numrank(Jval))
            end
        end

        m_tp = d_in * (d_in + 1) ÷ 2
        m_total = m_tp + n_Jz
        bp_bound = bp_rank_real(m_total)

        @printf("  %3.0f   | %3d  | %3d   | %5d | %5d  | %5d | %4d | %4d  %s\n",
                Jmax, d_in, d_det, m_tp, n_Jz, m_total, bp_bound,
                max_rank, max_rank <= bp_bound ? "✓" : "✗")
    end
end


# ============================================================================
# Summary: scattering complexity interpretation
# ============================================================================
function scattering_summary()
    println("\n" * "="^72)
    println("Physical interpretation: scattering complexity")
    println("="^72)
    println("""
The Barvinok-Pataki bound for inclusive scattering channels says:

  r_Kraus ≤ √(d_in² + k)

where:
  - d_in  = input Hilbert space dimension (incoming particles × spin × partial waves)
  - k     = number of conservation laws imposed as linear constraints
  - r     = Kraus rank of the inclusive channel (# of distinct scattering pathways)

Physical meaning:
  - Without conservation laws: r ≤ d_in (Choi's theorem). The channel can use
    all d_in pathways — no simplification.
  - Each conservation law (J_z, parity, charge, ...) adds constraints but also
    adds to m, potentially allowing slightly higher rank in principle.
  - However, symmetry typically REDUCES the effective dimension by block-
    diagonalising the Choi matrix, which in practice gives much lower ranks.
  - The BP bound provides a UNIVERSAL upper bound on scattering complexity
    that depends only on the number of constraints, not on the details
    of the interaction.
""")
end


# ============================================================================
# Run all examples
# ============================================================================

println("Barvinok-Pataki for Inclusive Scattering Channels")
println("=" ^ 52)
println("Solver: $solver_available")

example1_generic_scattering()
example2_electron_electron()
example3_compton()
scattering_summary()
