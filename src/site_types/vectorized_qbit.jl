# Additional Qubit gates
# ----------------------
ITensors.op(::OpName"Id", ::SiteType"Qubit") = [
    1 0
    0 1
]
ITensors.op(::OpName"id", st::SiteType"Qubit") = ITensors.op(OpName("Id"), st)

function ITensors.op(::OpName"CH", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("Id"), st)
    h = op(OpName("H"), st)
    return kron(proj0, id) + kron(proj1, h)
end

function ITensors.op(::OpName"CPhase", st::SiteType"Qubit"; ϕ::Real)
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("Id"), st)
    phase = op(OpName("Phase"), st; ϕ=ϕ)
    return kron(proj0, id) + kron(proj1, phase)
end

function ITensors.op(::OpName"CU3", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    u = ITensors.op(OpName("U"), st; θ=θ, ϕ=ϕ, λ=λ)
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("id"), st)
    return kron(proj0, id) + kron(proj1, u)
end

function ITensors.op(::OpName"CU", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    return ITensors.op(OpName("CU3"), st; θ=θ, ϕ=ϕ, λ=λ)
end

function ITensors.op(::OpName"CU1", st::SiteType"Qubit"; λ::Real)
    return ITensors.op(OpName("CU3"), st; θ=0, ϕ=0, λ=λ)
end

function ITensors.op(::OpName"CCCCNOT", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("Id"), st)
    not = op(OpName("X"), st)
    return kron(proj0, proj0, proj0, proj0, id) +
           kron(proj0, proj0, proj0, proj1, id) +
           kron(proj0, proj0, proj1, proj0, id) +
           kron(proj0, proj0, proj1, proj1, id) +
           kron(proj0, proj1, proj0, proj0, id) +
           kron(proj0, proj1, proj0, proj1, id) +
           kron(proj0, proj1, proj1, proj0, id) +
           kron(proj0, proj1, proj1, proj1, id) +
           kron(proj1, proj0, proj0, proj0, id) +
           kron(proj1, proj0, proj0, proj1, id) +
           kron(proj1, proj0, proj1, proj0, id) +
           kron(proj1, proj0, proj1, proj1, id) +
           kron(proj1, proj1, proj0, proj0, id) +
           kron(proj1, proj1, proj0, proj1, id) +
           kron(proj1, proj1, proj1, proj0, id) +
           kron(proj1, proj1, proj1, proj1, not)
end

"""
    ITensors.space(st::SiteType"vQubit"; dim = 2)

Create the Hilbert space for a site of type "vQubit", i.e. a vectorised
qbit, where the vectorisation is performed wrt the Pauli transfer matrix
basis of `Mat(ℂ²)`, composed of the identity matrix and the three Pauli matrices.
"""
ITensors.space(::SiteType"vQubit") = 4

# "Pauli transfer matrix" basis
# -----------------------------
function ptmbasis(nqbits::Int)
    st = SiteType("Qubit")
    σ =
        (1 / sqrt(2)) .*
        [Matrix(I, 2, 2), ITensors.op("X", st), ITensors.op("Y", st), ITensors.op("Z", st)]

    if nqbits == 1
        return σ  # and don't bother with the rest
    else
        # σxn = Base.product(s, s, s... [n times])
        σxn = Base.product(repeat([σ], nqbits)...)
        # Each element of σxn is a list (s_1, s_2, ..., s_n) where s_i is a Pauli matrix.
        # For example, for two bits, if `tp(a, b) = "$a ⊗ $b" then
        #   S = ["s$n" for n in 0:3];
        #   [tp(s...) for s in Base.product(repeat([S], 2)...)] -->
        #    "s0 ⊗ s0"  "s0 ⊗ s1"  "s0 ⊗ s2"  "s0 ⊗ s3"
        #    "s1 ⊗ s0"  "s1 ⊗ s1"  "s1 ⊗ s2"  "s1 ⊗ s3"
        #    "s2 ⊗ s0"  "s2 ⊗ s1"  "s2 ⊗ s2"  "s2 ⊗ s3"
        #    "s3 ⊗ s0"  "s3 ⊗ s1"  "s3 ⊗ s2"  "s3 ⊗ s3"

        # By calling `Base.vec` on this matrix we stack its columns, but we want to unroll
        # by rows instead, so we transpose it first. The call to `permutedims` below
        # does this transposition "on all dimensions".
        # FIXME Why do we need rows instead of columns here?
        tensorproducts = [kron(s...) for s in σxn]
        perm = reverse(ntuple(i -> i, Val{nqbits}()))  # (nqbits, nqbits - 1, ..., 2, 1)
        tensorproducts_transposed = permutedims(tensorproducts, perm)

        # We transform each of them in the tensor product s_1 ⊗ s_2 ⊗ ... ⊗ s_n.
        return Base.vec(tensorproducts_transposed)
    end
end

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vQubit")
    v = ITensors.state(sn, SiteType("Qubit"))
    return LindbladVectorizedTensors.vec(kron(v, v'), ptmbasis(1))
end
function vop(sn::StateName, ::SiteType"vQubit")
    sn = LindbladVectorizedTensors.statenamestring(sn)
    on = sn[1] == 'v' ? sn[2:end] : sn
    return LindbladVectorizedTensors.vec(
        LindbladVectorizedTensors.try_op(OpName(on), SiteType("Qubit")), ptmbasis(1)
    )
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"0", st::SiteType"vQubit") = vstate(sn, st)
ITensors.state(sn::StateName"1", st::SiteType"vQubit") = vstate(sn, st)

# States (vectorized operators)
# -----------------------------
ITensors.state(sn::StateName"vId", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"vX", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"vY", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"vZ", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"vH", st::SiteType"vQubit") = vop(sn, st)

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vQubit")
    d = Int(log2(size(mat, 1)))
    return LindbladVectorizedTensors.vec(x -> mat * x, ptmbasis(d))
end
function postmultiply(mat, ::SiteType"vQubit")
    d = Int(log2(size(mat, 1)))
    return LindbladVectorizedTensors.vec(x -> x * mat, ptmbasis(d))
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the S=1/2 site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vQubit"; kwargs...)
    name = strip(String(ITensors.name(on))) # Remove extra whitespace
    if name == "Id"
        return Matrix(1.0I, 4, 4)
    end
    dotloc = findfirst("⋅", name)
    # This returns the position of the cdot in the operator name String.
    # It is `nothing` if no cdot is found.
    if !isnothing(dotloc)
        on1, on2 = nothing, nothing
        on1 = name[1:prevind(name, dotloc.start)]
        on2 = name[nextind(name, dotloc.start):end]
        # If the OpName `on` is written correctly, i.e. it is either "A⋅" or "⋅A" for some
        # A, then either `on1` or `on2` has to be empty (not both, not neither of them).
        if (on1 == "" && on2 == "") || (on1 != "" && on2 != "")
            throw(
                ArgumentError(
                    "Invalid operator name: $name. Operator name is not \"Id\" or of the " *
                    "form \"A⋅\" or \"⋅A\"",
                ),
            )
        end
        # name == "⋅A" -> on1 is an empty string
        # name == "A⋅" -> on2 is an empty string
        if on1 == ""
            mat = LindbladVectorizedTensors.try_op(
                OpName(on2), SiteType("Qubit"); kwargs...
            )
            return postmultiply(mat, st)
        elseif on2 == ""
            mat = LindbladVectorizedTensors.try_op(
                OpName(on1), SiteType("Qubit"); kwargs...
            )
            return premultiply(mat, st)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
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

"""
    adjointmap(X::OpName; kwargs...)

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
