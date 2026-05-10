# Space of spin-1/2 particles (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vS=1/2"; dim = 2)

Create the Hilbert space for a site of type "vS=1/2", i.e. a vectorised
spin-1/2 particle, where the vectorisation is performed wrt the generalised
Gell-Mann basis of `Mat(ℂ²)`, composed of Hermitian traceless matrices
together with the identity matrix.
"""
ITensors.space(::SiteType"vS=1/2") = 4

# Elements of and operators on Mat(ℂ²) are expanded wrt the basis {Λᵢ}ᵢ₌₁⁴ of
# generalised Gell-Mann matrices (plus a multiple of the identity).
# An element A ∈ Mat(ℂ²) is representeb by the a vector v such that
#     vᵢ = tr(Λᵢ A),
# while a linear map L : Mat(ℂ²) → Mat(ℂ²) by the matrix ℓ such that
#     ℓᵢⱼ = tr(Λᵢ L(Λⱼ)).

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vS=1/2")
    v = ITensors.state(sn, SiteType("S=1/2"))
    return _hilbertschmidt_vec(kron(v, v'), gellmannbasis(2))
end
function vop(sn::StateName, ::SiteType"vS=1/2")
    return _hilbertschmidt_vec(op(statenamestring(sn), siteind("S=1/2")), gellmannbasis(2))
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"Up", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"Dn", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"↑", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"↓", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"X+", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"X-", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"Y+", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"Y-", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"Z+", st::SiteType"vS=1/2") = vstate(sn, st)
ITensors.state(sn::StateName"Z-", st::SiteType"vS=1/2") = vstate(sn, st)

# States representing vectorised operators
# ----------------------------------------
ITensors.state(sn::StateName"Sx", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"Sy", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"Sz", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"X", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"Y", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"Z", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"σx", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"σy", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"σz", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"Id", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"N", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"S+", st::SiteType"vS=1/2") = vop(sn, st)
ITensors.state(sn::StateName"S-", st::SiteType"vS=1/2") = vop(sn, st)
