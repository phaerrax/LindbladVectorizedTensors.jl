"""
    statenamestring(sn::StateName{T}) where T

Return the name of the state `sn` as a string.
"""
function statenamestring(sn::StateName{T}) where {T}
    return string(T)
end

"""
    sitenumber(i::Index)

Return the site number of the given Index, i.e. the number in the "n=N" tag.
"""
function sitenumber(i::Index)
    st = split(replace(string.(tags.(i)), "\"" => ""), ",")
    sitetag = first(filter!(s -> occursin("n=", s), st))
    # TODO what if sitetag is empty?
    siten = replace(sitetag, "n=" => "")
    return parse(Int, siten)
end

"""
    jwstring(; start, stop, op::AbstractString="F")

Return the vector ``(op, start + 1, op, start + 2, ..., op, stop - 1)``.
"""
function jwstring(; start, stop, op::AbstractString="F")
    return collect(Iterators.flatten([(op, start + k) for k in 1:(stop - start - 1)]))
end

# Vectorisation utilities
# =======================
# In order to use tr(x'*y) as a tool to extract coefficient the basis must
# of course be orthonormal wrt this inner product.
# The canonical basis or the Gell-Mann one are okay.

"""
    vec(A::Matrix, basis::Vector)

Compute the vector of coefficients of the matrix `A` wrt the basis `basis`.
"""
function vec(A::Matrix, basis::Vector)
    return [tr(b' * A) for b in basis]
end

"""
    vec(L::Function, basis::Vector)

Compute the matrix of coefficients of the linear map `L` wrt the basis `basis`.
The linearity of the map is not checked, so using this function with non-linear
functions leads to undefined results.
"""
function vec(L::Function, basis::Vector)
    return [tr(bi' * L(bj)) for (bi, bj) in Base.product(basis, basis)]
end

"""
    partialtrace(sites::Vector{Index{Int64}}, v::MPS, j::Int)

Compute the partial trace, on the `j`-th site of `sites`, of the matrix
represented, as vectorised, by the MPS `v`.
The result is a `Vector` containing the coordinates of the partial trace
(i.e. the result is still in vectorised form).
"""
function partialtrace(sites::Vector{Index{Int64}}, v::MPS, j::Int)
    # Custom contraction sequence: we start from v[end] and we contract it with
    # a vecId state; we contract the result with v[end-1] and so on, until we
    # get to v[j].  # Then we do the same starting from v[1].
    x = ITensor(1.0)
    for i in length(sites):-1:(j + 1)
        x = v[i] * state(sites[i], "vecId") * x
    end
    y = ITensor(1.0)
    for i in 1:(j - 1)
        y = v[i] * state(sites[i], "vecId") * y
    end
    z = y * v[j] * x
    # Now `vector(z)` is the coordinate vector of the partial trace.
    return vector(z)
end

# Generalised Gell-Mann matrices
# ==============================
"""
    gellmannmatrix(j, k, dim)

Return the `(j,k)` generalised Gell-Mann matrix of dimension `dim`, normalised
wrt the Hilbert-Schmidt inner product ``(A,B) = tr(A†B)``.
The matrices are indexed as follows:

    * if ``j > k`` the matrix is symmetric and traceless;
    * if ``j < k`` the matrix is antisymmetric;
    * if ``j = k`` the matrix is diagonal.

In particular, ``j = k = dim`` gives a matrix proportional to the identity.
The two indices `j` and `k` determine the non-zero coefficients of the matrix.
The whole set of (different) Gell-Mann matrices that can be generated with this
function is a basis of ``Mat(ℂᵈⁱᵐ)``.
"""
function gellmannmatrix(j, k, dim)
    if j > dim || k > dim || j < 0 || k < 0
        throw(DomainError)
    end
    m = zeros(ComplexF64, dim, dim)
    if j > k
        m[j, k] = 1 / sqrt(2)
        m[k, j] = 1 / sqrt(2)
    elseif k > j
        m[j, k] = -im / sqrt(2)
        m[k, j] = im / sqrt(2)
    elseif j == k && j < dim
        for i in 1:j
            m[i, i] = 1
        end
        m[j + 1, j + 1] = -j
        m .*= sqrt(1 / (j * (j + 1)))
    else
        for i in 1:dim
            m[i, i] = 1 / sqrt(dim)
        end
    end
    return m
end

"""
    gellmannbasis(dim)

Return a list containing the "Hermitian basis" of ``Mat(ℂᵈⁱᵐ)``, i.e. composed
of the ``dim²`` generalised Gell-Mann matrices.
"""
function gellmannbasis(dim)
    return [gellmannmatrix(j, k, dim) for (j, k) in [Base.product(1:dim, 1:dim)...]]
    # We need to splat the result from `product` so that the result is a list
    # of matrices (a Vector) and not a Matrix.
end

"""
    canonicalmatrix(i, j, dim)

Return the (`i`,`j`) element of the canonical basis of ``Mat(ℂᵈⁱᵐ)``, i.e. a
`dim`×`dim` matrix whose element on the `i`-th row and `j`-th column is ``1``,
and zero elsewhere.
"""
function canonicalmatrix(i, j, dim)
    m = zeros(ComplexF64, dim, dim)
    m[i, j] = 1
    return m
end

"""
    canonicalbasis(dim)

Return a list of the matrices in the canonical basis of ``Mat(ℂᵈⁱᵐ)``. 
The list is ordered corresponding to column-based vectorisation, i.e.

    canonicalbasis(dim)[j] = canonicalmatrix((j-1)%dim + 1, (j-1)÷dim + 1, dim)

with ``j ∈ {1,…,dim²}``. With this ordering,
``vec(A)ⱼ = tr(canonicalbasis(dim)[j]' * A)``.
"""
function canonicalbasis(dim)
    return [canonicalmatrix(i, j, dim) for (i, j) in [Base.product(1:dim, 1:dim)...]]
end

# (Von Neumann) Entropy
# =====================
"""
    vonneumannentropy(ψ::MPS, sites::Vector{Index{Int64}}, n::Int)

Compute the entanglement entropy of the biparition ``(1,…,n)|(n+1,…,N)`` of
the system in state described by the MPS `ψ` (defined on the sites `sites`),
using its Schmidt decomposition.
"""
function vonneumannentropy(ψ₀::MPS, sites::Vector{Index{Int64}}, n::Int)
    ψ = orthogonalize(ψ₀, n)
    # Decompose ψ[n] in singular values, treating the Link between sites n-1 and n
    # and the physical index as "row index"; the remaining index, the Link
    # between n and n+1, is the "column index".
    _, S, _ = svd(ψ[n], (linkind(ψ, n - 1), sites[n]))
    # Compute the square of the singular values (aka the Schmidt coefficients
    # of the bipartition), and from them the entropy.
    sqdiagS = [S[j, j]^2 for j in ITensors.dim(S, 1)]
    return -sum(p -> p * log(p), sqdiagS; init=0.0)
end

# Chop (from Mathematica)
# =======================
# Imitation of Mathematica's "Chop" function

"""
    chop(x::Real; tolerance=1e-10)

Truncates `x` to zero if it is less than `tolerance`.
"""
function chop(x::Real; tolerance=1e-10)
    return abs(x) > tolerance ? x : zero(x)
end

"""
    chop(x::Complex; tolerance=1e-10)

Truncates the real and/or the imaginary part of `x` to zero if they are less
than `tolerance`.
"""
function chop(x::Complex; tolerance=1e-10)
    return Complex(chop(real(x)), chop(imag(x)))
end

# MPS and MPO utilities
# =====================
"""
    chain(left::MPS, right::MPS)

Concatenate `left` and `right`, returning `left` ``⊗`` `right`.
"""
function chain(left::MPS, right::MPS)
    # This function is like mpnum's `chain`: it takes two MPSs ans concatenates
    # them. The code is "inspired" from `ITensors.MPS` in mps.jl:308.

    midN = length(left) # The site with the missing link between the two MPSs.
    # First of all we shift the Link tags of the given MPSs, so that the final
    # enumeration of the tags is correct.
    # Note that in each MPS the numbers in the Link tags do not follow the
    # numbering of the Sites on which it is based, they always start from 1.
    for j in eachindex(left)
        replacetags!(left[j], "l=$j", "l=$j"; tags="Link")
        replacetags!(left[j], "l=$(j-1)", "l=$(j-1)"; tags="Link")
    end
    for j in eachindex(right)
        replacetags!(right[j], "l=$j", "l=$(midN+j)"; tags="Link")
        replacetags!(right[j], "l=$(j-1)", "l=$(midN+j-1)"; tags="Link")
    end
    # "Shallow" concatenation of the MPSs (there's still a missing link).
    M = MPS([left..., right...])
    # We create a "trivial" Index of dimension 1 and add it to the two sites
    # which are not yet connected.
    # The Index has dimension 1 because this is a tensor product between the
    # states so there's no correlation between them.
    missing_link = Index(1; tags="Link,l=$midN")
    M[midN] = M[midN] * state(missing_link, 1)
    M[midN + 1] = state(dag(missing_link), 1) * M[midN + 1]

    return M
end

"""
    chain(left::MPO, right::MPO)

Concatenate `left` and `right`, returning `left` ``⊗`` `right`.
"""
function chain(left::MPO, right::MPO)
    # Like the previous `chain`, but for MPOs.
    midN = length(left)# The site with the missing link between the two MPOs.
    for j in eachindex(right)
        replacetags!(right[j], "l=$j", "l=$(midN+j)"; tags="Link")
        replacetags!(right[j], "l=$(j-1)", "l=$(midN+j-1)"; tags="Link")
    end
    M = MPO([left..., right...])
    missing_link = Index(1; tags="Link,l=$midN")
    M[midN] = M[midN] * state(missing_link, 1)
    M[midN + 1] = M[midN + 1] * state(dag(missing_link), 1)
    # The order of the Indexes in M[midN] and M[midN+1] ns not what we would get
    # if we built an MPO on the whole system as usual; namely, the two "Link"
    # Indexes are swapped.
    # This however should not matter since ITensor does not care about the order.
    return M
end

# Varargs versions
"""
    chain(a::MPS, b...)

Concatenate the given MPSs into a longer MPS, returning their tensor product.
"""
chain(a::MPS, b...) = chain(a, chain(b...))

"""
    chain(a::MPO, b...)

Concatenate the given MPOs into a longer MPO, returning their tensor product.
"""
chain(a::MPO, b...) = chain(a, chain(b...))

"""
    embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPO)

Embed `slice`, defined on a subset `range` of `sites`, into a MPO which covers
the whole `sites`.

The MPO is extended by filling the empty spots with an "Id" operator, therefore
an operator with OpName "Id" is required to be defined for the SiteTypes of the
remaining sites.

# Arguments
- `sites::Array{Index{Int64}}`: the sites of the whole system.
- `range::UnitRange{Int}`: the range spanned by `slice`.
- `slice::MPO`: the MPO to be extended.
"""
function embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPO)
    # TODO: compute automatically on which sites the MPO is defined, without
    # having to supply the range explicitly as an argument.
    if length(slice) != length(range)
        throw(DimensionMismatch("slice and range must have the same size."))
    end
    if !issubset(range, eachindex(sites))
        throw(BoundsError(range, sites))
    end

    if range[begin] == 1 && range[end] == length(sites)
        mpo = slice
    elseif range[begin] == 1
        mpo = chain(slice, MPO(sites[(range[end] + 1):end], "Id"))
    elseif range[end] == length(sites)
        mpo = chain(MPO(sites[1:(range[begin] - 1)], "Id"), slice)
    else
        mpo = chain(
            MPO(sites[1:(range[begin] - 1)], "Id"),
            slice,
            MPO(sites[(range[end] + 1):end], "Id"),
        )
    end
    return mpo
end

"""
    embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPS)

Embed `slice`, defined on a subset `range` of `sites`, into a MPS which covers
the whole `sites` (to be interpreted as a vectorised operator).

The MPS is extended by filling the empty spots with a "vecId" operator,
therefore an operator with OpName "vecId" is required to be defined for the
SiteTypes of the remaining sites.

# Arguments
- `sites::Array{Index{Int64}}`: the sites of the whole system.
- `range::UnitRange{Int}`: the range spanned by `slice`.
- `slice::MPS`: the MPS to be extended.
"""
function embed_slice(sites::Array{Index{Int64}}, range::UnitRange{Int}, slice::MPS)
    if length(slice) != length(range)
        throw(DimensionMismatch("slice e range must have the same size."))
    end
    if !issubset(range, eachindex(sites))
        throw(BoundsError(range, sites))
    end

    if range[begin] == 1 && range[end] == length(sites)
        mpo = slice
    elseif range[begin] == 1
        mpo = chain(slice, MPS(sites[(range[end] + 1):end], "vecId"))
    elseif range[end] == length(sites)
        mpo = chain(MPS(sites[1:(range[begin] - 1)], "vecId"), slice)
    else
        mpo = chain(
            MPS(sites[1:(range[begin] - 1)], "vecId"),
            slice,
            MPS(sites[(range[end] + 1):end], "vecId"),
        )
    end
    return mpo
end

# Other utilities
# ===============
"""
    consecutivepairs(v::AbstractVector)

Return a list of Strings "(a,b)" formed by all adjacent items in `v`.
"""
function consecutivepairs(v::AbstractVector)
    return string.("(", v[1:(end - 1)], ",", v[2:end], ")")
end

"""
    try_op(on::OpName, st::SiteType; kwargs...)

Return the matrix of an ITensor operator, if it exists, trying first the
```julia
op(::OpName, ::SiteType; kwargs...)
```
syntax, and then
```julia
op!(::ITensor, ::OpName, ::SiteType, ::Index...; kwargs...)
```
if the former returns nothing.
"""
function try_op(on::OpName, st::SiteType; kwargs...)
    stname = String(ITensors.tag(st))
    opstring = String(ITensors.name(on))
    # Try calling a function of the form:
    #    op(::OpName, ::SiteType; kwargs...)
    # which returns a Julia matrix
    mat = ITensors.op(on, st; kwargs...)
    if isnothing(mat)
        # Otherwise try calling a function of the form
        #    op!(::ITensor, ::OpName, ::SiteType, ::Index...; kwargs...)
        dummy = siteind(stname)
        Op = ITensor(prime(dummy), ITensors.dag(dummy))
        r = ITensors.op!(Op, on, st, dummy; kwargs...)
        if isnothing(r)
            throw(
                ArgumentError(
                    "Overload of \"op\" or \"op!\" functions not found for " *
                    "operator name \"$opstring\" and Index tag $(tags(dummy)).",
                ),
            )
        end
        mat = matrix(Op)
    end
    if isnothing(mat)
        throw(
            ArgumentError(
                "Overload of \"op\" or \"op!\" functions not found for operator " *
                "name \"$opstring\" and Index tag \"$stname\".",
            ),
        )
    end
    return mat
end

"""
    try_op(on::OpName, st::SiteType, d::Int; kwargs...)

Like try_op(on::OpName, st::SiteType; kwargs...) (see [`try_op`](@ref)), but with an
additional Int argument so that it can be used by SiteTypes without a fixed dimension.
"""
function try_op(on::OpName, st::SiteType, d::Int; kwargs...)
    stname = String(ITensors.tag(st))
    opstring = String(ITensors.name(on))
    # Try calling a function of the form:
    #    op(::OpName, ::SiteType, ::Int; kwargs...)
    # which returns a Julia matrix
    mat = ITensors.op(on, st, d; kwargs...)
    if isnothing(mat)
        # Otherwise try calling a function of the form
        #    op!(::ITensor, ::OpName, ::SiteType, ::Index..., ::Int; kwargs...)
        dummy = siteind(stname)
        Op = ITensor(prime(dummy), ITensors.dag(dummy))
        r = ITensors.op!(Op, on, st, dummy, d; kwargs...)
        if isnothing(r)
            throw(
                ArgumentError(
                    "Overload of \"op\" or \"op!\" functions not found for " *
                    "operator name \"$opstring\" and Index tag $(tags(dummy)).",
                ),
            )
        end
        mat = matrix(Op)
    end
    if isnothing(mat)
        throw(
            ArgumentError(
                "Overload of \"op\" or \"op!\" functions not found for operator " *
                "name \"$opstring\" and Index tag \"$stname\".",
            ),
        )
    end
    return mat
end
