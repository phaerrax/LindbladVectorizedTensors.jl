using LindbladVectorizedTensors, ITensors, ITensorMPS, LinearAlgebra
using Test

function trace(x::MPS)
    vid = MPS(siteinds(x), "vId")
    return dot(vid, x)
end

function gmat(v; dim=2)
    sum(vi * b for (vi, b) in zip(v, LindbladVectorizedTensors.gellmannbasis(dim)))
end

@testset "Definition of S=1/2 states" begin
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

# Test that trace(..., x_vec) == dot(..., x_vec)
# Sx_exp_vec = [trace(apply(op("Sy⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]
# Sy_exp_vec = [trace(apply(op("Sy⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]
# Sz_exp_vec = [trace(apply(op("Sy⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]
# Test that the vectorization is done correctly by measuring the expectation values of some
# observables on a random state.

function expect_vec(x::MPS, name::AbstractString)
    return [dot(MPS(siteinds(x), j -> j == n ? name : "Id"), x) for n in 1:length(x)]
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

    @testset "Qubit" begin
        sites = siteinds("Qubit", 4)
        x = random_mps(sites; linkdims=4)
        x_vec = vec_projector(x)
        sites_vec = siteinds(x_vec)

        y_exp = expect(x, "Y")
        y_exp_vec = [trace(apply(op("Y⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]

        h_exp = expect(x, "H")
        h_exp_vec = [trace(apply(op("H⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]

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
