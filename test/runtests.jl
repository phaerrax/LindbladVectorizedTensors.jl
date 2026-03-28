using LindbladVectorizedTensors, ITensors, ITensorMPS, LinearAlgebra
using Test

include("vectorization.jl")

function gmat(v; dim=2)
    sum(vi * b for (vi, b) in zip(v, LindbladVectorizedTensors.gellmannbasis(dim)))
end

@testset "vS=1/2 states" begin
    vs = siteind("vS=1/2")
    s = siteind("S=1/2")
    @test state(vs, "↑") == state(vs, "Up")
    @test state(vs, "↓") == state(vs, "Dn")
    @test gmat(vector(state(vs, "X+"))) ≈ 1/2 * (I + matrix(op(s, "σx")))
    @test gmat(vector(state(vs, "X-"))) ≈ 1/2 * (I - matrix(op(s, "σx")))
    @test gmat(vector(state(vs, "Y+"))) ≈ 1/2 * (I + matrix(op(s, "σy")))
    @test gmat(vector(state(vs, "Y-"))) ≈ 1/2 * (I - matrix(op(s, "σy")))
end

@testset "Qubit MPO -> vQubit MPS vectorization" begin
    @test all(vectorization_expvalues(; atol=1e-14))
end
