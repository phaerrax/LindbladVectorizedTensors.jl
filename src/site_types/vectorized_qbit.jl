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
    return _hilbertschmidt_vec(kron(v, v'), ptmbasis(1))
end
function vop(sn::StateName, ::SiteType"vQubit")
    return _hilbertschmidt_vec(op(statenamestring(sn), siteind("Qubit")), ptmbasis(1))
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
