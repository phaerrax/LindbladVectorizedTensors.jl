using LindbladVectorizedTensors, ITensors, ITensorMPS, LinearAlgebra
using Test

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
    sum(vi * b for (vi, b) in zip(v, LindbladVectorizedTensors.gellmannbasis(dim)))
end

@testset "Left- and right-multiplication operators" begin
    N = 4
    sites = siteinds("vS=1/2", N)
    x = random_mps(sites; linkdims=4)
    # Note that we need to generate _real_ MPSs, since the states are Hermitian matrices
    # hence linear combinations of Gell-Mann matrices with real coefficients.
    # The tests fail if the MPS is complex.

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
        @test scalar(state(vs, "Id") * apply(op("N⋅", vs), ρT)) ≈ 1/expm1(β * ω)
    end
end

@testset "Vectorisation of operators" verbose=true begin
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
    end

    @testset "Qubit" begin
        sites = siteinds("Qubit", 4)
        x = random_mps(sites; linkdims=4)
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

        @test h_exp ≈ h_exp_vec
        @test h_exp ≈ h_exp_vec
        @test cp_exp ≈ cp_exp_vec
        @test ccx_exp ≈ ccx_exp_vec
    end
end
