export adjointmap_itensor

"""
    adjointmapmatrix(X::OpName; kwargs...)

Return the matrix (in the Gell-Mann basis) representing the action of ``X`` on a
state ``ρ`` as ``XρX⁻¹``. The argument must be a valid ITensors OpName for the Qubit site
type.
Additional parameters needed to specify the operator may be passed as keyword arguments.
"""
function adjointmapmatrix(X::OpName; kwargs...)
    op = ITensors.op(X, SiteType("Qubit"); kwargs...)
    # This is some square matrix. We can deduce the dimension, hence the number of qbits
    # on which X acts, from it.
    n_qbits = Int(log2(size(op, 1)))
    return LindbladVectorizedTensors.vec(a -> op * a * op', ptmbasis(n_qbits))
end

"""
    adjointmap_itensor(X::OpName, s1::Index, s_tail::Index...; kwargs...)

Return the ITensor representing the action of ``X`` on a state ``ρ`` as ``XρX⁻¹``, where
`X` acts on sites given by the Index list.
The argument must be a valid ITensors OpName for the Qubit site type.
Additional parameters needed to specify the operator may be passed as keyword arguments.
"""
function adjointmap_itensor(on::OpName, s1::Index, s_tail::Index...; kwargs...)
    # Adapted from the ITensors.op function for the `vOsc` site type in LindbladVectorizedTensors
    rs = reverse((s1, s_tail...))
    opmat = adjointmapmatrix(on; kwargs...)
    return ITensors.itensor(opmat, prime.(rs)..., dag.(rs)...)
end

function adjointmap_itensor(x::AbstractString, s1::Index, s_tail::Index...; kwargs...)
    return adjointmap_itensor(OpName(x), s1, s_tail...; kwargs...)
end

"""
    adjointmap_itensor(on::OpName, sites::Vector{<:Index}, n::Int...; kwargs...)

Return the ITensor representing the action of ``X`` on a state ``ρ`` as ``XρX⁻¹``, where
`X` acts on the given sites.
The argument must be a valid ITensors OpName for the Qubit site type.
Additional parameters needed to specify the operator may be passed as keyword arguments.
"""
function adjointmap_itensor(on::OpName, sites::Vector{<:Index}, n::Int...; kwargs...)
    s1, s_tail... = [sites[j] for j in n]
    return adjointmap_itensor(on, s1, s_tail...; kwargs...)
end

function adjointmap_itensor(x::AbstractString, sites::Vector{<:Index}, n::Int...; kwargs...)
    return adjointmap_itensor(OpName(x), sites, n...; kwargs...)
end

"""
    adjointmap_itensor(sites::Vector{<:Index}, on::OpName, n::Int...; kwargs...)

Return the ITensor representing the action of ``X`` on a state ``ρ`` as ``XρX⁻¹``, where
`X` acts on the given sites.
The argument must be a valid ITensors OpName for the Qubit site type.
Additional parameters needed to specify the operator may be passed as keyword arguments.
"""
function adjointmap_itensor(sites::Vector{<:Index}, on::OpName, n::Int...; kwargs...)
    return adjointmap_itensor(on, sites, n...; kwargs...)
end

function adjointmap_itensor(sites::Vector{<:Index}, x::AbstractString, n::Int...; kwargs...)
    return adjointmap_itensor(sites, OpName(x), n...; kwargs...)
end
