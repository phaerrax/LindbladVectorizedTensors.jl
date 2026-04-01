using ITensors.SiteTypes: _sitetypes

"""
    space(::SiteType"vFDot3")

Create the Hilbert space for a site of type "vFDot3", i.e. a mixed state describing a
site with a fermionic three-level quantum dot.
The density matrix is represented in the generalised Gell-Mann basis, composed
of Hermitian traceless matrices together with the identity matrix.

No conserved symmetries and quantum number labels are provided for this space.
"""
function ITensors.space(::SiteType"vFDot3")
    return (2^3)^2
end

# An element A ∈ Mat(ℂ⁴) is representeb by the a vector v such that
#     vᵢ = tr(Λᵢ A),
# while a linear map L : Mat(ℂ⁴) → Mat(ℂ⁴) by the matrix ℓ such that
#     ℓᵢⱼ = tr(Λᵢ L(Λⱼ)).

# Shorthand notation:
function vstate(sn::AbstractString, ::SiteType"vFDot3")
    v = ITensors.state(StateName(sn), SiteType("FDot3"))
    return LindbladVectorizedTensors.vec(kron(v, v'), gellmannbasis(2^3))
end
function vop(on::AbstractString, ::SiteType"vFDot3")
    return LindbladVectorizedTensors.vec(
        try_op(OpName(statenamestring(sn)), SiteType("FDot3")), gellmannbasis(2^3)
    )
end

# basis order:
# e_1 -> |∅⟩
# e_2 -> c₁†|∅⟩
# e_3 -> c₂†|∅⟩
# e_4 -> c₁† c₂†|∅⟩
# e_5 -> c₃†|∅⟩
# e_6 -> c₁† c₃†|∅⟩
# e_7 -> c₂† c₃†|∅⟩
# e_8 -> c₁† c₂† c₃†|∅⟩
ITensors.state(sn::StateName"Emp", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(::StateName"0", st::SiteType"vFDot3") = state(StateName("Emp"), st)
ITensors.state(::StateName"Vac", st::SiteType"vFDot3") = state(StateName("Emp"), st)
ITensors.state(::StateName"Vacuum", st::SiteType"vFDot3") = state(StateName("Emp"), st)

ITensors.state(sn::StateName"1", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(sn::StateName"2", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(sn::StateName"12", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(sn::StateName"3", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(sn::StateName"13", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(sn::StateName"23", st::SiteType"vFDot3") = vstate(sn, st)
ITensors.state(sn::StateName"123", st::SiteType"vFDot3") = vstate(sn, st)

ITensors.state(sn::StateName"Id", st::SiteType"vFDot3") = vop(sn, st)
ITensors.state(sn::StateName"n1", st::SiteType"vFDot3") = vop(sn, st)
ITensors.state(sn::StateName"n2", st::SiteType"vFDot3") = vop(sn, st)
ITensors.state(sn::StateName"n3", st::SiteType"vFDot3") = vop(sn, st)
ITensors.state(sn::StateName"ntot", st::SiteType"vFDot3") = vop(sn, st)

# Operator dispatch
# =================
function premultiply(mat, ::SiteType"vFDot3")
    return LindbladVectorizedTensors.vec(x -> mat * x, gellmannbasis(2^3))
end
function postmultiply(mat, ::SiteType"vFDot3")
    return LindbladVectorizedTensors.vec(x -> x * mat, gellmannbasis(2^3))
end

# The goal here is to define operators "A⋅" and "⋅A" in an automatic way whenever the
# OpName "A" is defined for the FDot3 site type.
# This is handy, but unless we find a better way to define this function this means that
# _every_ operator has to be written this way; we cannot just return op(on, st) at the end
# if no "⋅" is found, otherwise an infinite loop would be entered.
# We make an exception, though, for "Id" since it is an essential operator, and something
# would probably break if it weren't defined.
function ITensors.op(on::OpName, st::SiteType"vFDot3"; kwargs...)
    name = strip(String(ITensors.name(on))) # Remove extra whitespace
    if name == "Id"
        return Matrix(1.0I, 2^6, 2^6)
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
            mat = try_op(OpName(on2), SiteType("FDot3"); kwargs...)
            return postmultiply(mat, st)
        elseif on2 == ""
            mat = try_op(OpName(on1), SiteType("FDot3"); kwargs...)
            return premultiply(mat, st)
        else
            # This should logically never happen but, just in case, we throw an error.
            error("Unknown error with operator name $name")
        end
    else
        error("Operator name $name is not \"Id\" or of the form \"A⋅\" or \"⋅A\"")
    end
end

function dot_hamiltonian(::SiteType"vFDot3", energies, coulomb_repulsion, sitenumber::Int)
    E = OpSum()
    for k in 1:3
        E += energies[k] * gkslcommutator("n$k", sitenumber)
    end

    N² = gkslcommutator("ntot^2", sitenumber)
    N = gkslcommutator("ntot", sitenumber)

    return E + 0.5coulomb_repulsion * (N² - N)
end

function exchange_interaction(::SiteType"vFDot3", s1::Index, s2::Index)
    stypes1 = _sitetypes(s1)
    stypes2 = _sitetypes(s2)
    if (SiteType("vFDot3") in stypes1) && !(SiteType("vFDot3") in stypes2)
        return exchange_interaction(st, sitenumber(s1), sitenumber(s2))
    elseif (SiteType("vFDot3") in stypes2) && !(SiteType("vFDot3") in stypes1)
        return exchange_interaction(st, sitenumber(s2), sitenumber(s1))
    else
        # Return an error if no implementation is found for any type.
        throw(
            ArgumentError("No vFDot3 site type found in either $(tags(s1)) or $(tags(s2)).")
        )
    end
end

function exchange_interaction(::SiteType"vFDot3", dot_site::Int, other_site::Int)
    return sum([
        gkslcommutator("c†$k", dot_site, "σ-", other_site) +
        gkslcommutator("c$k", dot_site, "σ+", other_site) for k in 1:3
    ])
end
