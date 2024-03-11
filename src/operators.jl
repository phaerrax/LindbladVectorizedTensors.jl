# Easy construction of the unitary part of the GKSL equation
# ----------------------------------------------------------
function makeopsumpairs(args...)
    isodd(length(args)) &&
        throw(error("Cannot make pairs out of an odd number of arguments."))
    return [(args[i], args[i + 1]) for i in eachindex(args)[1:2:(end - 1)]]
end

"""
    gkslcommutator(operators::Tuple{String,Int}...)

Given one or more tuples `(a1, n1), (a2, n2), …`, return an OpSum representing the operation
``x ↦ -i[A_1 A_2 …, x]`` where ``A_i`` is an operator which consists of `ai` at site `ni`
and the identity elsewhere. The string `a1` must be an existing ITensor OpName whose
variants `a1⋅` and `⋅a1` are defined (this function, however, doesn't perform any checks).

# Examples
```julia-repl
julia> gkslcommutator(("σ+", 1), ("σ-", 2))
sum(
  0.0 - 1.0im σ+⋅(1,) σ-⋅(2,)
  0.0 + 1.0im ⋅σ+(1,) ⋅σ-(2,)
)
```
"""
function gkslcommutator(operators::Tuple{String,Int}...)
    l = OpSum()
    opnames = first.(operators)
    sites = last.(operators)
    l += (-im, collect(Iterators.flatten(zip(opnames .* "⋅", sites)))...)
    l += (+im, collect(Iterators.flatten(zip("⋅" .* opnames, sites)))...)
    return l
end

"""
    gkslcommutator(args...)

Given a sequence `a1, n1, a2, n2, …` where each `ai` is a `String` and each `ni` is an
`Int`, return an OpSum representing the operation ``x ↦ -i[A_1 A_2 …, x]`` where ``A_i``
is an operator which consists of `ai` at site `ni` and the identity elsewhere. The string
`a1` must be an existing ITensor OpName whose variants `a1⋅` and `⋅a1` are defined
(this function, however, doesn't perform any checks).

# Examples
```julia-repl
julia> gkslcommutator("σ+", 1, "σ-", 2)
sum(
  0.0 - 1.0im σ+⋅(1,) σ-⋅(2,)
  0.0 + 1.0im ⋅σ+(1,) ⋅σ-(2,)
)
```
"""
gkslcommutator(args...) = gkslcommutator(makeopsumpairs(args...)...)

"""
    gkslcommutator_itensor(sites::Vector{<:Index}, operators::Tuple{String,Int}...)

Given a vector of ITensors site indices `sites` and a sequence `a1, n1, a2, n2, …` where
each `ai` is a `String` and each `ni` is an `Int`, return an ITensor with the site indices
`s[n1]`, `s[n1]'`, `s[n2]`, `s[n2]'` and so on which represents the operation
``x ↦ -i[A_1 A_2 …, x]`` where ``A_i`` is an operator consisting of `ai` at site
`ni` and the identity elsewhere. Each `ai` string must be an existing ITensors OpName whose
variants `a1⋅` and `⋅a1` are defined (this function, however, doesn't perform any checks).


# Examples

We start from a system of three 1/2-spins:
```julia
sites = siteinds("vS=1/2", 3)
```

Take an operator ``U`` which is ``Sˣ`` on the second site and the identity on the
others. The commutator ``ρ ↦ -i[U,ρ]`` is given by
```julia
gkslcommutator_itensor(sites, "Sx", 2)
```

An operator ``V`` which is ``Sy`` on the first site, the identity on the second one, and
``Sz`` on the third one. The commutator ``ρ ↦ -i[V,ρ]`` is given by
```julia
gkslcommutator_itensor(sites, "Sy", 1, "Sz", 3)
```
"""
function gkslcommutator_itensor(sites::Vector{<:Index}, operators::Tuple{String,Int}...)
    operator_names = first.(operators)
    site_numbers = last.(operators)
    if !allunique(sites)
        error("Some sites are repeated in the list. Please use unique site indices.")
        # This possibility is not allowed for now, since it would create issues in the
        # multiplication loop below. Basically, if two ITensors with the same indices
        # are multiplied together (with a simple `*`) then both indices get contracted
        # and we end up with a scalar. We should use `apply` instead, but it does not
        # work with OneITensor objects...
    end
    lmult = ITensors.OneITensor()
    rmult = ITensors.OneITensor()
    for (on, j) in zip(operator_names, site_numbers)
        lmult *= op("$on⋅", sites, j)
        rmult *= op("⋅$on", sites, j)
    end
    return -im * (lmult - rmult)
end

function gkslcommutator_itensor(sites::Vector{<:Index}, args...)
    return gkslcommutator_itensor(sites, makeopsumpairs(args...)...)
end
