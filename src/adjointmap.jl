export adjointmap_itensor

"""
    ilog(base::Int, arg::Int)

Integer logarithm: the largest integer `n` such that `base^n ≤ arg`.

julia> isqrt(4, 18)
2

julia> isqrt(3, 81)
4

"""
function ilog(base::Int, arg::Int)
    if base ≤ 1
        throw(DomainError("logarithm base must be greater than 1"))
    end
    if arg < 1
        throw(DomainError("logarithm argument must be greater or equal than 1"))
    end
    l = 0
    x = arg
    while x ≥ base
        l += 1
        x = x/base
    end
    return l
end

"""
    adjointmap_matrix(op::AbstractMatrix; dim)

Return the matrix (in the Gell-Mann basis) representing the action of an operator ``X`` on a
state ``ρ`` as ``XρX⁻¹``, given the matrix representation of ``X``.
"""
function adjointmap_matrix(op::AbstractMatrix; dim, site_type)
    # `op` is some square matrix. We can deduce the dimension, hence the number of qbits
    # on which it acts, from its size.
    nsites = ilog(dim, size(op, 1))
    vbasis = vectorizationbasis(site_type, nsites; dim)
    return _hilbertschmidt_vec(a -> op * a * op', vbasis)
end

"""
    adjointmap_itensor(t::ITensor, orig_sites::Vector{<:Index}, vec_sites::Vector{<:Index})
    adjointmap_itensor(op_name, sites::Index...; kwargs...)
    adjointmap_itensor(op_name, sites::Vector{<:Index}, n::Int...; kwargs...)

Return the ITensor representing the action of an operator ``X`` on a state ``ρ`` as
``XρX⁻¹``, where ``X`` acts on the given site(s).

``X`` may be given as an ITensor already acting on the (unvectorized) sites `orig_sites`, or
as a valid ITensors operator name `op_name` for the Qubit site type, in which case
additional parameters needed to specify the operator may be passed as keyword arguments.
"""
function adjointmap_itensor(
    t::ITensor, orig_sites::Vector{<:Index}, vec_sites::Vector{<:Index}
)
    # qubit_sites and vqubit_sites must be given in the SAME site order, i.e.,
    # qubit_sites[i] and vqubit_sites[i] must refer to the same physical site --- the
    # combiner does not care what tags/ids those indices carry beyond their dimension.
    cmb = combiner(orig_sites...)
    cmb_index = combinedind(cmb)
    mat = matrix(cmb' * t * cmb, cmb_index', cmb_index)
    # matrix(T, i, j) returns a matrix M such that T = itensor(M, i, j), that is,
    #   M[a, b] = T[i => a, j => b]
    # so with this line we are making sure that cmb_index' indexes rows and cmb_index
    # indexes columns.

    # Check if the original sites are supported by this package, and retrieve the dimension
    # of the relative Hilbert space.
    # The dimension must be computed from the indices themselves, since it's not uniquely
    # determined by the site type (e.g. Bosons, whose dimension is determined by the user).
    # At the same time we need the site type name so that we can call the appropriate
    # vectorisation function.
    stypes = ITensors.SiteTypes._sitetypes.(orig_sites)  # all tags
    common_stypes = intersect(stypes...)  # keep only shared tags
    filter!(!=(SiteType("Site")), common_stypes)  # remove "Site", which is usually there
    common_stype = only(common_stypes)  # only one must remain

    dims = dim.(orig_sites)
    if !allequal(dims)
        error("adjointmap_itensor only supports operators defined on a single site type")
    end

    op_mat = adjointmap_matrix(mat; dim=first(dims), site_type=common_stype)
    return ITensors.itensor(op_mat, prime.(vec_sites)..., dag.(vec_sites)...)
end

function adjointmap_itensor(
    on::Union{AbstractString,OpName}, s1::Index, s_tail::Index...; kwargs...
)
    vqubit_sites = [s1, s_tail...]
    # We need to create a temporary list of "Qubit" sites to use with the ITensor-based
    # adjointmap_itensor method.
    qubit_sites = [Index(2, "Qubit") for _ in vqubit_sites]
    t = ITensors.op(on, qubit_sites...; kwargs...)
    return adjointmap_itensor(t, qubit_sites, vqubit_sites)
end

# Variant with site list and indices given separately.
function adjointmap_itensor(
    on::Union{AbstractString,OpName}, sites::Vector{<:Index}, n::Int...; kwargs...
)
    s1, s_tail... = [sites[j] for j in n]
    return adjointmap_itensor(on, s1, s_tail...; kwargs...)
end

# Variant with the list of sites given before the operator name. We don't really expose this
# syntax, but ITensor has it, so maybe there's a point in defining it.
function adjointmap_itensor(
    sites::Vector{<:Index}, on::Union{AbstractString,OpName}, n::Int...; kwargs...
)
    return adjointmap_itensor(on, sites, n...; kwargs...)
end
