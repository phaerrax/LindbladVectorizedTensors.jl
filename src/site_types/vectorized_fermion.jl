# Space of spin-1/2 particles (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vFermion")

Create the Hilbert space for a site of type "vFermion", i.e. a vectorised
spin-1/2 particle, where the vectorisation is performed wrt the generalised
Gell-Mann basis of `Mat(ℂ²)`, composed of Hermitian traceless matrices
together with the identity matrix.
"""
ITensors.space(::SiteType"vFermion") = 4

# Elements of and operators on Mat(ℂ²) are expanded wrt the basis {Λᵢ}ᵢ₌₁⁴ of
# generalised Gell-Mann matrices (plus a multiple of the identity).
# An element A ∈ Mat(ℂ²) is representeb by the a vector v such that
#     vᵢ = tr(Λᵢ A),
# while a linear map L : Mat(ℂ²) → Mat(ℂ²) by the matrix ℓ such that
#     ℓᵢⱼ = tr(Λᵢ L(Λⱼ)).

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vFermion")
    v = ITensors.state(sn, SiteType("Fermion"))
    return _hilbertschmidt_vec(kron(v, v'), gellmannbasis(2))
end
function vop(sn::StateName, ::SiteType"vFermion")
    return _hilbertschmidt_vec(
        op(statenamestring(sn), siteind("Fermion")), gellmannbasis(2)
    )
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"Emp", st::SiteType"vFermion") = vstate(sn, st)
ITensors.state(sn::StateName"Occ", st::SiteType"vFermion") = vstate(sn, st)

function ITensors.state(::StateName"Up", st::SiteType"vFermion")
    return ITensors.state(StateName("Occ"), st)
end
function ITensors.state(::StateName"Dn", st::SiteType"vFermion")
    return ITensors.state(StateName("Emp"), st)
end

# States representing vectorised operators
# ----------------------------------------
ITensors.state(sn::StateName"Id", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"N", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"F", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"A", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"a", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"Adag", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"adag", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"A†", st::SiteType"vFermion") = vop(sn, st)
ITensors.state(sn::StateName"a†", st::SiteType"vFermion") = vop(sn, st)
