export gkslcommutator, gkslcommutator_itensor

# Easy construction of the unitary part of the GKSL equation
# ----------------------------------------------------------
function makeopsumpairs(args...)
    isodd(length(args)) &&
        throw(error("Cannot make pairs out of an odd number of arguments."))
    return [(args[i], args[i + 1]) for i in eachindex(args)[1:2:(end - 1)]]
end

"""
    gkslcommutator(ops::Tuple{AbstractString,Int}...)
    gkslcommutator(ops...)

Given one or more tuples `(a1, n1), (a2, n2), …`, or alternatively a sequence `a1, n1, a2,
n2, …` where each `ai` is a `String` and each `ni` is an `Int`, return an OpSum representing
the operation ``x ↦ -i[A_1 A_2 …, x]`` where ``A_i`` is an operator which consists of `ai`
on the site `ni` and the identity elsewhere. Each string `ai` must be the name of an
existing ITensor operator (this function however does not perform any check).

# Examples
```julia-repl
julia> gkslcommutator(("σ+", 1), ("σ-", 2))
sum(
  0.0 - 1.0im σ+⋅(1,) σ-⋅(2,)
  0.0 + 1.0im ⋅σ+(1,) ⋅σ-(2,)
)

julia> gkslcommutator("A", 3, "Adag", 5)
sum(
  0.0 - 1.0im A⋅(3,) Adag⋅(5,)
  0.0 + 1.0im ⋅A(3,) ⋅Adag(5,)
)

```
"""
function gkslcommutator end

function gkslcommutator(operators::Tuple{AbstractString,Int}...)
    l = OpSum()
    opnames = first.(operators)
    sites = last.(operators)
    l += (-im, collect(Iterators.flatten(zip(opnames .* "⋅", sites)))...)
    l += (+im, collect(Iterators.flatten(zip("⋅" .* opnames, sites)))...)
    return l
end

gkslcommutator(args...) = gkslcommutator(makeopsumpairs(args...)...)

"""
    gkslcommutator_itensor(s::Vector{<:Index}, ops::Tuple{String,Int}...)
    gkslcommutator_itensor(s::Vector{<:Index}, ops...)

Given a vector of ITensors site indices `s` and a sequence `a1, n1, a2, n2, …`, where each
`ai` is a string and each `ni` is an integer, return an ITensor with site indices
`s[n1]`, `s[n1]'`, `s[n2]`, `s[n2]'` and so on, which represents the operation ``x ↦ -i[A_1
A_2 …, x]``, where ``A_i`` is an operator consisting of `ai` acting on site `ni` and the
identity elsewhere. Each `ai` string must be the name of an existing ITensor operator of
the respective non-vectorised site type(s).

# Examples

We start from a system of three 1/2-spins:

```julia-repl
julia> sites = siteinds("vS=1/2", 3)
3-element Vector{Index{Int64}}:
 (dim=4|id=...|"Site,n=1,vS=1/2")
 (dim=4|id=...|"Site,n=2,vS=1/2")
 (dim=4|id=...|"Site,n=3,vS=1/2")

```

Take an operator ``U`` which is made up of ``Sx`` on the second site and the identity on the
others. The commutator ``ρ ↦ -i[U,ρ]`` is given by

```julia-repl
julia> gkslcommutator_itensor(sites, "Sx", 2)
ITensor ord=2 (dim=4|id=...|"Site,n=2,vS=1/2")' (dim=4|id=...|"Site,n=2,vS=1/2")
NDTensors.Dense{ComplexF64, Vector{ComplexF64}}

```

An operator ``V`` which is ``Sy`` on the first site, the identity on the second one, and
``Sz`` on the third one. The commutator ``ρ ↦ -i[V,ρ]`` is given by

```julia-repl
julia> gkslcommutator_itensor(sites, "Sy", 1, "Sz", 3)
ITensor ord=4 (dim=4|id=...|"Site,n=1,vS=1/2")' (dim=4|id=...|"Site,n=1,vS=1/2") (dim=4|id=...|"Site,n=3,vS=1/2")' (dim=4|id=...|"Site,n=3,vS=1/2")
NDTensors.Dense{ComplexF64, Vector{ComplexF64}}

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

#=
FIXME This functions does not work unless the tuple with the keyword arguments is empty (in
which case this method is kind of useless): define for example

julia> ITensors.op(::OpName"rX", st::SiteType"Qubit"; x) = cis(x) * op(OpName("X"), st)

and try running:

julia> gkslcommutator_itensor(vs, ("rX", (x=4)), 3)
ERROR: syntax: invalid named tuple element ""rX"" around REPL[20]:1
Stacktrace:
 [1] top-level scope
   @ REPL[20]:1

We need to rethink how to handle keyword arguments. For now let's disable this variant.

===============
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
=#

function gkslcommutator_itensor(sites::Vector{<:Index}, args...)
    return gkslcommutator_itensor(sites, makeopsumpairs(args...)...)
end
