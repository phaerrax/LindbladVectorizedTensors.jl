export gkslcommutator, gkslcommutator_itensor

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

#=
Some gates require parameters. We need to pass these parameters to gkslcommutator_itensor
as well, if we want to build the commutator.

Example 1: a commutator with the U gate
---------------------------------------
We need to get to the expression (with some random values for the angles)
    -im * (op("U⋅", s[n]; θ=pi/4, ϕ=0, λ=pi/8) - op("⋅U", s[n]; θ=pi/4, ϕ=0, λ=pi/8)).
The parameters are best associated to the gate name, rather than to the site, therefore
it makes sense to redefine gkslcommutator_itensor so that it accepts arguments as
    gkslcommutator_itensor(sites, X, n, X', n'...)
where X, X',... contain information about the gate name _and_ its parameters.
Since the parameters are given to ITensors.op as keyword arguments, i.e. a NamedTuple,
it might make sense to specify X as a Tuple with the operator name (a String) in the first
position and the parameters (a NamedTuple) in the second position:
    X = ("U", (θ=pi/4, ϕ=0, λ=pi/8))
and later build the relevant operators using
    on, args = X
    lmult = ITensors.op("$on⋅", s, n; args...)
    rmult = ITensors.op("⋅$on", s, n; args...)
    -im * (lmult - rmult)

Example 2: a commutator with the H gate
---------------------------------------
What if the gate doesn't have any parameter? Well, luckily the above syntax works even if
`args` is empty, but we still need to provide it, so we will use
    X = ("H", ())
and then go on as above; `args` will be an empty tuple, and we have
    ITensors.op("H⋅", s, n'; ()...) == ITensors.op("H⋅", s, n').

Example 3: both the above gates at the same time
------------------------------------------------
Note that the `gkslcommutator_itensor(::Vector{<:Index}, ::Tuple{String,Int}...)` function
(the original one) works like this: you pass the vector of site indices and then some
of {String,Int} pairs, which get slurped into a Vector of Pair{String,Int} items:
    gkslcommutator_itensor(s, ("A", 1), ("B", 3)).
For the purposes of this library, "A" and "B" need to be replaced by the X and X' things
above, so by (another!) tuple:
    gkslcommutator_itensor(s, (("U", (θ=pi/4, ϕ=0, λ=pi/8)), 1), (("H", ()), 3)).
The second and third argument of `gkslcommutator_itensor` here (try this on the REPL!) are
Tuple{Tuple{String, Any}, Int} items. The new signature should maybe then be
    gkslcommutator_itensor(::Vector{<:Index}, ::Tuple{Tuple{String, Any},Int}...).
=#

function gkslcommutator_itensor(
    sites::Vector{<:Index}, items::Tuple{Tuple{String,Any},Int}...
)
    # Unpacking the arguments:
    # - separate operator names/args and site numbers
    operators = first.(items)
    site_numbers = last.(items)
    if !allunique(site_numbers)
        error("Some sites are repeated in the list. Please use unique site indices.")
        # This possibility is not allowed for now, since it would create issues in the
        # multiplication loop below. Basically, if two ITensors with the same indices
        # are multiplied together (with a simple `*`) then both indices get contracted
        # and we end up with a scalar. We should use `apply` instead, but it does not
        # work with OneITensor objects...
    end
    # - separate operator names from parameters
    operator_names = first.(operators)
    operator_args = last.(operators)

    lmult = ITensors.OneITensor()
    rmult = ITensors.OneITensor()
    for (on, args, j) in zip(operator_names, operator_args, site_numbers)
        lmult *= op("$on⋅", sites, j; args...)
        rmult *= op("⋅$on", sites, j; args...)
    end
    return -im * (lmult - rmult)
end

#function gkslcommutator_itensor(sites::Vector{<:Index}, args...)
#    return gkslcommutator_itensor(sites, makeopsumpairs(args...)...)
#end
