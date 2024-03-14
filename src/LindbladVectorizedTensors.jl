module LindbladVectorizedTensors

using ITensors
using IterTools
using LinearAlgebra

import ITensors.state, ITensors.op, ITensors.space

export chain,
    chop,
    consecutivepairs,
    embed_slice,
    gellmannbasis,
    gellmannmatrix,
    jwstring,
    sitenumber,
    vec,
    vonneumannentropy

include("utils.jl")

include("site_types/spinhalf.jl")
include("site_types/vectorized_spinhalf.jl")
include("site_types/fermion.jl")
include("site_types/vectorized_fermion.jl")
include("site_types/oscillator.jl")
include("site_types/electron.jl")
include("site_types/vectorized_electron.jl")
include("site_types/fermionic_dot3.jl")
include("site_types/vectorized_fermionic_dot3.jl")

export dissipator_loss, dissipator_gain, dissipator, mixedlindbladplus, mixedlindbladminus

include("site_types/vectorized_oscillator.jl")

export gkslcommutator, gkslcommutator_itensor

include("operators.jl")

include("spin_chain.jl")
include("exchange_interaction.jl")

end
