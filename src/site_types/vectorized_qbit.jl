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
    return LindbladVectorizedTensors.vec(
        try_op(OpName(statenamestring(sn)), SiteType("Qubit")), ptmbasis(1)
    )
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"0", st::SiteType"vQubit") = vstate(sn, st)
ITensors.state(sn::StateName"1", st::SiteType"vQubit") = vstate(sn, st)

# States (vectorized operators)
# -----------------------------
ITensors.state(sn::StateName"Id", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"X", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"Y", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"Z", st::SiteType"vQubit") = vop(sn, st)
ITensors.state(sn::StateName"H", st::SiteType"vQubit") = vop(sn, st)

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
