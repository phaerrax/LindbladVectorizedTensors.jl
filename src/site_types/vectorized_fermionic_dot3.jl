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

# basis order:
# e_1 -> |∅⟩
# e_2 -> c₁†|∅⟩
# e_3 -> c₂†|∅⟩
# e_4 -> c₁† c₂†|∅⟩
# e_5 -> c₃†|∅⟩
# e_6 -> c₁† c₃†|∅⟩
# e_7 -> c₂† c₃†|∅⟩
# e_8 -> c₁† c₂† c₃†|∅⟩

ITensors.state(::StateName"0", st::SiteType"vFDot3") = state(StateName("Emp"), st)
ITensors.state(::StateName"Vac", st::SiteType"vFDot3") = state(StateName("Emp"), st)
ITensors.state(::StateName"Vacuum", st::SiteType"vFDot3") = state(StateName("Emp"), st)

# States derived from the FDot3 site type
register_vectorized_names(
    SiteType("vFDot3");
    states=("Emp", "1", "2", "12", "3", "13", "23", "123"),
    operators=("Id", "n1", "n2", "n3", "ntot"),
)

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
