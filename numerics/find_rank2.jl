#!/usr/bin/env julia
#
# find_rank2.jl — Find cases where the SDP optimizer genuinely has rank > 1
#
# Strategy: pick constraint values in the INTERIOR of the joint numerical
# range (not achievable by any pure state).  Then every feasible state
# must be mixed, and the extreme-point optimizer has rank ≥ 2.
#

include("hydrogen_sdp.jl")

using LinearAlgebra

"""
    pure_state_range(basis, operators)

Compute the joint numerical range of `operators` over pure states on the shell.
Returns a matrix where row i is (⟨ψᵢ|A₁|ψᵢ⟩, …, ⟨ψᵢ|Aₖ|ψᵢ⟩) for random pure |ψᵢ⟩.
"""
function pure_state_range(basis, operators; num_samples=10000)
    d = length(basis)
    k = length(operators)
    samples = zeros(num_samples, k)
    for i in 1:num_samples
        # Random pure state (Haar measure)
        ψ = randn(ComplexF64, d)
        ψ ./= norm(ψ)
        for (j, A) in enumerate(operators)
            samples[i, j] = real(dot(ψ, A * ψ))
        end
    end
    return samples
end

function main()
    println("=" ^ 70)
    println("Searching for rank-2 optimizers")
    println("=" ^ 70)

    # Work on n=2 shell (dim=8) for clarity
    n = 2
    basis = shell_basis(n)
    d = length(basis)
    Jz = Jz_op(basis)
    J2 = Matrix(J2_op(basis))
    L2 = Matrix(L2_op(basis))
    H_fs = fine_structure_op(basis)

    println("\nShell n=$n, dim=$d")
    println("\nStep 1: Map pure-state joint numerical range for (Jz, J², L²)")

    ops = [Jz, J2, L2]
    samples = pure_state_range(basis, ops, num_samples=50000)

    println("  Jz range:  [$(minimum(samples[:,1])), $(maximum(samples[:,1]))]")
    println("  J² range:  [$(minimum(samples[:,2])), $(maximum(samples[:,2]))]")
    println("  L² range:  [$(minimum(samples[:,3])), $(maximum(samples[:,3]))]")

    # Find interior point: average of pure-state samples
    # (centroid of numerical range is always interior)
    centroid = vec(mean(samples, dims=1))
    println("\n  Centroid of joint numerical range:")
    println("    ⟨Jz⟩ = $(centroid[1])")
    println("    ⟨J²⟩ = $(centroid[2])")
    println("    ⟨L²⟩ = $(centroid[3])")

    # Check: is the centroid achievable by a pure state?
    # Find closest pure state to centroid
    dists = [norm(samples[i,:] - centroid) for i in 1:size(samples,1)]
    min_dist = minimum(dists)
    println("    Closest pure state distance: $min_dist")

    println("\nStep 2: Solve SDPs with interior constraint values")
    println("=" ^ 70)

    # Try the centroid as constraint values
    B = 1.0

    # m = 3 constraints (BP bound: r ≤ 2)
    println("\n--- m=3: Jz + J² + L² at centroid values ---")
    cons3 = [("Jz", Jz, centroid[1]),
             ("J²", J2, centroid[2]),
             ("L²", L2, centroid[3])]
    E3, r3, eigs3, bp3 = solve_hydrogen_sdp(n, B, cons3)
    println("  Energy: $E3")
    println("  Rank:   $r3  (BP bound: $bp3)")
    println("  Eigs:   $(round.(eigs3[1:min(5,d)], digits=6))")

    # Try with Jz = 0 (symmetric point)
    println("\n--- m=3: Jz=0, J²=centroid, L²=centroid ---")
    cons3b = [("Jz", Jz, 0.0),
              ("J²", J2, centroid[2]),
              ("L²", L2, centroid[3])]
    E3b, r3b, eigs3b, bp3b = solve_hydrogen_sdp(n, B, cons3b)
    println("  Energy: $E3b")
    println("  Rank:   $r3b  (BP bound: $bp3b)")
    println("  Eigs:   $(round.(eigs3b[1:min(5,d)], digits=6))")

    # m=2: Jz=0 + J²=centroid (BP bound: r ≤ 1)
    println("\n--- m=2: Jz=0, J²=centroid ---")
    cons2 = [("Jz", Jz, 0.0),
             ("J²", J2, centroid[2])]
    E2, r2, eigs2, bp2 = solve_hydrogen_sdp(n, B, cons2)
    println("  Energy: $E2")
    println("  Rank:   $r2  (BP bound: $bp2)")
    println("  Eigs:   $(round.(eigs2[1:min(5,d)], digits=6))")

    # Now try with SPECIFIC interior values designed to force rank > 1
    # Key insight: if we constrain J² to a value between the eigenvalues
    # (3/4 and 15/4), AND constrain Jz=0, AND constrain L² to a non-integer
    # value, the pure-state constraints become overdetermined
    println("\n--- Designed interior points ---")
    for (jz_val, j2_val, l2_val) in [
        (0.0, 2.0, 1.0),   # J² between 3/4 and 15/4
        (0.0, 1.5, 0.5),
        (0.0, 2.5, 1.5),
        (0.0, 3.0, 1.0),
        (0.5, 2.0, 1.0),
    ]
        cons = [("Jz", Jz, jz_val),
                ("J²", J2, j2_val),
                ("L²", L2, l2_val)]
        E, r, eigs, bp = solve_hydrogen_sdp(n, B, cons)
        status = r > 1 ? "★ RANK>1" : ""
        println(@sprintf("  Jz=%.1f, J²=%.1f, L²=%.1f → rank=%d (BP≤%d) E=%.4f %s",
            jz_val, j2_val, l2_val, r, bp, E, status))
        if r > 1
            println("    Eigenvalues: $(round.(eigs[1:min(5,d)], digits=6))")
        end
    end

    # Try on n=3 shell for more room
    println("\n\n--- n=3 shell (dim=18), designed interior points ---")
    n3 = 3
    basis3 = shell_basis(n3)
    d3 = length(basis3)
    Jz3 = Jz_op(basis3)
    J2_3 = Matrix(J2_op(basis3))
    L2_3 = Matrix(L2_op(basis3))

    ops3 = [Jz3, J2_3, L2_3]
    samples3 = pure_state_range(basis3, ops3, num_samples=50000)
    centroid3 = vec(mean(samples3, dims=1))
    println("  Centroid: Jz=$(round(centroid3[1],digits=3)), J²=$(round(centroid3[2],digits=3)), L²=$(round(centroid3[3],digits=3))")

    for (jz_val, j2_val, l2_val) in [
        (centroid3[1], centroid3[2], centroid3[3]),
        (0.0, 3.0, 2.0),
        (0.0, 4.0, 3.0),
        (0.0, 2.0, 1.5),
        (0.0, 5.0, 3.0),
    ]
        cons = [("Jz", Jz3, jz_val),
                ("J²", J2_3, j2_val),
                ("L²", L2_3, l2_val)]
        E, r, eigs, bp = solve_hydrogen_sdp(n3, B, cons)
        status = r > 1 ? "★ RANK>1" : ""
        println(@sprintf("  Jz=%.2f, J²=%.2f, L²=%.2f → rank=%d (BP≤%d) E=%.4f %s",
            jz_val, j2_val, l2_val, r, bp, E, status))
        if r > 1
            println("    Eigenvalues: $(round.(eigs[1:min(5,d3)], digits=6))")
        end
    end

    println("\n\nDone.")
end

using Statistics
main()
