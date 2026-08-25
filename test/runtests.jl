using Test, Documenter, LindbladVectorizedTensors
using ITensors, ITensorMPS, LinearAlgebra

@testset "Documentation examples" begin
    doctest(LindbladVectorizedTensors; manual=false)
end

trace(x::MPS) = dot(MPS(siteinds(x), "Id"), x)

function expect_trace(x::MPS, name::AbstractString)
    # Valore atteso di A sullo stato ρ calcolato con ⟨V(1), V(A⋅) V(ρ)⟩
    return [trace(apply(op(name * "⋅", siteinds(x), n), x)) for n in 1:length(x)]
end

function expect_vec(x::MPS, name::AbstractString)
    # Valore atteso di A sullo stato ρ calcolato con ⟨V(ρ), V(A)⟩
    return [
        dot(x, MPS(ComplexF64, siteinds(x), j -> j == n ? name : "Id")) for n in 1:length(x)
    ]
end

function gmat(v; dim=2)
    # Reconstruct the matrix by its coefficients (in the Gell-Mann basis).
    return sum(vi * b for (vi, b) in zip(v, LindbladVectorizedTensors.gellmannbasis(dim)))
end

@testset "Left- and right-multiplication operators" begin
    # For an orthonormal basis {b_i} (wrt the Hilbert-Schmidt inner product) the
    # coefficient vector of a matrix M is v_i = tr(b_i' * M).
    # `premultiply(t, vst)` must be the matrix representation, in that same basis, of the
    # linear map ρ ↦ t*ρ, and likewise `postmultiply(t, vst)` for ρ ↦ ρ*t.
    coeffs(m, basis) = [tr(b' * m) for b in basis]

    function test_multiply(st_name; dim=nothing)
        st = SiteType(st_name)
        vst = SiteType("v" * st_name)
        site_idx, site_dim = if isnothing(dim)
            siteinds(st_name, 1), ITensors.space(st)
        else
            siteinds(st_name, 1; dim=dim), dim
        end
        basis = LindbladVectorizedTensors.vectorizationbasis(st, 1; dim=site_dim)

        # We don't really need to generate physically meaningful operators and states here
        # to check that pre/postmultiply work: random matrices will do.

        t_mat = rand(ComplexF64, site_dim, site_dim)
        t = ITensors.op(t_mat, site_idx)

        ρ_mat = rand(ComplexF64, site_dim, site_dim)
        ρ_vec = coeffs(ρ_mat, basis)

        return LindbladVectorizedTensors.premultiply(t, vst) * ρ_vec ≈
               coeffs(t_mat * ρ_mat, basis) &&
               LindbladVectorizedTensors.postmultiply(t, vst) * ρ_vec ≈
               coeffs(ρ_mat * t_mat, basis)
    end

    @test test_multiply("S=1/2")
    @test test_multiply("Electron")
    @test test_multiply("Fermion")
    @test test_multiply("Qubit")
    # Here we check the case where the dimension must be explicitly provided.
    @test test_multiply("Boson"; dim=3) && test_multiply("Boson"; dim=4)

    @testset "MPS interoperability" begin
        N = 4
        sites = siteinds("vS=1/2", N)
        x = random_mps(sites; linkdims=4)
        # Note that we need to generate _real_ MPSs, since the states are Hermitian matrices
        # hence linear combinations of Gell-Mann matrices with real coefficients.
        # The tests might fail if the MPS is complex (but why though? TODO find out!).

        @test expect_trace(x, "Sx") ≈ expect_vec(x, "Sx")
        @test expect_trace(x, "Sy") ≈ expect_vec(x, "Sy")
        @test expect_trace(x, "Sz") ≈ expect_vec(x, "Sz")

        sites = siteinds("vBoson", N; dim=5)
        x = random_mps(sites; linkdims=4)

        @test expect_trace(x, "N") ≈ expect_vec(x, "N")
        @test expect_trace(x, "X") ≈ expect_vec(x, "X")
        @test expect_trace(x, "A") ≈ expect_vec(x, "A")

        sites = siteinds("vFermion", N)
        x = random_mps(sites; linkdims=4)

        @test expect_trace(x, "N") ≈ expect_vec(x, "N")
        @test expect_trace(x, "A") ≈ expect_vec(x, "A")
    end
end

@testset "Definition of vectorised states" verbose=true begin
    @testset "S=1/2" begin
        vs = siteind("vS=1/2")
        s = siteind("S=1/2")
        @test state(vs, "↑") == state(vs, "Up")
        @test state(vs, "↓") == state(vs, "Dn")

        @test gmat(vector(state(vs, "X+"))) ≈ 1/2 * (I + matrix(op(s, "σx")))
        @test gmat(vector(state(vs, "X-"))) ≈ 1/2 * (I - matrix(op(s, "σx")))
        @test gmat(vector(state(vs, "Y+"))) ≈ 1/2 * (I + matrix(op(s, "σy")))
        @test gmat(vector(state(vs, "Y-"))) ≈ 1/2 * (I - matrix(op(s, "σy")))
        @test gmat(vector(state(vs, "Z+"))) ≈ 1/2 * (I + matrix(op(s, "σz")))
        @test gmat(vector(state(vs, "Z-"))) ≈ 1/2 * (I - matrix(op(s, "σz")))
        @test gmat(vector(state(vs, "Z+"))) ≈ matrix(op(s, "ProjUp"))
        @test gmat(vector(state(vs, "Z-"))) ≈ matrix(op(s, "ProjDn"))
    end

    @testset "Boson" begin
        d = 5
        vs = siteind("vBoson"; dim=d)

        ρ = [state(vs, string(n-1)) for n in 1:d]
        for n in 1:(d - 1)
            @test apply(op("Adag⋅ * ⋅A", vs), ρ[n]) ≈ n * ρ[n + 1]
        end

        ω = 1/2 + rand()
        β = 1 + 10rand()
        ρT = state(vs, "ThermEq"; frequency=ω, temperature=1/β)
        # Watch out! The average number of bosons here is not exactly 1/expm1(β*ω) because
        # the Hilbert space is truncated at the d-th level.
        avgn = sum(n * exp(-β*ω*n) for n in 0:(d - 1)) / sum(exp(-β*ω*n) for n in 0:(d - 1))
        @test scalar(state(vs, "Id") * apply(op("N⋅", vs), ρT)) ≈ avgn
    end
end

@testset "Vectorisation of operators" verbose=true begin
    @testset "vec_projector keyword arguments" begin
        N = 2
        sites = siteinds("S=1/2", N)
        x = random_mps(ComplexF64, sites; linkdims=4)

        sites_vec = siteinds("vS=1/2", N)
        @test siteinds(vec_projector(x; existing_sites=sites_vec)) == sites_vec
        @test maxlinkdim(vec_projector(x; maxdim=3)) == 3
        @test maxlinkdim(vec_projector(x; cutoff=1e-4)) ≤ 4
    end

    @testset "S=1/2" begin
        N = 4
        sites = siteinds("S=1/2", N)
        x = random_mps(ComplexF64, sites; linkdims=4)
        x_vec = vec_projector(x)
        sites_vec = siteinds(x_vec)

        @test expect(x, "Sx") ≈ expect_vec(x_vec, "Sx")
        @test expect(x, "Sy") ≈ expect_vec(x_vec, "Sy")
        @test expect(x, "Sz") ≈ expect_vec(x_vec, "Sz")
    end

    @testset "Boson" begin
        N = 4
        sites = siteinds("Boson", N; dim=6)
        x = random_mps(ComplexF64, sites; linkdims=4)
        x_vec = vec_projector(x)
        sites_vec = siteinds(x_vec)

        @test expect(x, "N") ≈ expect_vec(x_vec, "N")
        @test expect(x, "X") ≈ expect_vec(x_vec, "X")
        @test expect(x, "A") ≈ expect_vec(x_vec, "A")

        adagb_exp = [
            dot(x, apply(op("a†b", sites, i, i + 1), x)) for
            i in 1:length(x) if i + 1 <= length(x)
        ]
        adagb_exp_vec = [
            trace(apply(op("ab†⋅", sites_vec, i, i + 1), x_vec)) for
            i in 1:length(x_vec) if i + 1 <= length(x_vec)
        ]
        # In `adagb_exp_vec` we must use the adjoint of the operator in `adagb_exp` for the
        # test to pass. This is the same issue we encounter in the "Left- and
        # right-multiplication operators" test set above.
        @test adagb_exp ≈ adagb_exp_vec
    end

    @testset "Qubit" begin
        sites = siteinds("Qubit", 4)
        x = random_mps(ComplexF64, sites; linkdims=4)
        x_vec = vec_projector(x)
        sites_vec = siteinds(x_vec)

        y_exp = expect(x, "Y")
        y_exp_vec = expect_trace(x_vec, "Y")

        h_exp = expect(x, "H")
        h_exp_vec = expect_trace(x_vec, "H")

        angle = pi * rand()
        cp_exp = [
            dot(x, apply(op("CPhase", sites, i, i + 2; ϕ=angle), x)) for
            i in 1:length(x) if i + 2 <= length(x)
        ]
        cp_exp_vec = [
            trace(apply(op("CPhase⋅", sites_vec, i, i + 2; ϕ=angle), x_vec)) for
            i in 1:length(x_vec) if i + 2 <= length(x_vec)
        ]

        ccx_exp = [
            dot(x, apply(op("CCNOT", sites, i, i + 1, i + 2), x)) for
            i in 1:length(x) if i + 2 <= length(x)
        ]
        ccx_exp_vec = [
            trace(apply(op("CCNOT⋅", sites_vec, i, i + 1, i + 2), x_vec)) for
            i in 1:length(x_vec) if i + 2 <= length(x_vec)
        ]

        @test y_exp ≈ y_exp_vec
        @test h_exp ≈ h_exp_vec
        @test cp_exp ≈ cp_exp_vec
        @test ccx_exp ≈ ccx_exp_vec
    end
end

@testset "Adjoint-map tensor" verbose=true begin
    sites = siteinds("Qubit", 3)
    v = random_mps(ComplexF64, sites; linkdims=4)
    v_vec = vec_projector(v)
    sites_vec = siteinds(v_vec)

    # Normal tensor. Check that it works with one index and with more than one.
    u = ITensors.op("RandomUnitary", sites[1])
    u_vec = adjointmap_itensor(u, [sites[1]], [sites_vec[1]])
    @test vec_projector(apply(u, v); existing_sites=sites_vec) ≈ apply(u_vec, v_vec)
    u = ITensors.op("RandomUnitary", sites[2], sites[3])
    u_vec = adjointmap_itensor(u, [sites[2], sites[3]], [sites_vec[2], sites_vec[3]])
    @test vec_projector(apply(u, v); existing_sites=sites_vec) ≈ apply(u_vec, v_vec)

    # Tensor with keyword arguments.
    θ = rand()
    u = ITensors.op("Ry", sites, 1; θ)
    u_vec = adjointmap_itensor("Ry", sites_vec, 1; θ)
    @test vec_projector(apply(u, v); existing_sites=sites_vec) ≈ apply(u_vec, v_vec)

    # Vector{<:Index}-based call forms.
    t_a = adjointmap_itensor("CX", sites_vec, 1, 2)
    t_b = adjointmap_itensor(sites_vec, "CX", 1, 2)
    t_c = adjointmap_itensor(
        ITensors.op("CX", sites, 1, 2), [sites[1], sites[2]], [sites_vec[1], sites_vec[2]]
    )
    @test t_a ≈ t_b ≈ t_c

    # Other site types than Qubit.
    sites = siteinds("S=1/2", 2)
    v = random_mps(ComplexF64, sites; linkdims=4)
    v_vec = vec_projector(v)
    sites_vec = siteinds(v_vec)
    xy = ITensors.op("X", sites[1]) * ITensors.op("Y", sites[2])
    xy_vec = adjointmap_itensor(xy, [sites[1], sites[2]], [sites_vec[1], sites_vec[2]])
    @test vec_projector(apply(xy, v); existing_sites=sites_vec) ≈ apply(xy_vec, v_vec)

    sites = siteinds("Boson", 3; dim=4)
    v = random_mps(ComplexF64, sites; linkdims=4)
    v_vec = vec_projector(v)
    sites_vec = siteinds(v_vec)
    aadag = ITensors.op("A", sites[1]) * ITensors.op("Adag", sites[2])
    aadag_vec = adjointmap_itensor(
        aadag, [sites[1], sites[2]], [sites_vec[1], sites_vec[2]]
    )
    # Note that `aadag` isn't unitary so we must pass `normalize=false` to `vec_projector`,
    # otherwise it will renormalise `apply(aadag, v)` before computing the projector.
    @test vec_projector(apply(aadag, v); existing_sites=sites_vec, normalize=false) ≈
        apply(aadag_vec, v_vec)
end
