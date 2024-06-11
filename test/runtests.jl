using LindbladVectorizedTensors, ITensors
using Test

include("vectorization.jl")

@testset "Qubit MPO -> vQubit MPS vectorization" begin
    @test all(vectorization_expvalues(; atol=1e-14))
end
