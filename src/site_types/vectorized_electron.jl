# Space of electrons (vectorised)
# ========================================
"""
    ITensors.space(st::SiteType"vElectron")

Create the Hilbert space for a site of type "vElectron", i.e. a mixed state describing a
site with a 1/2-spin degree of freedom.
The density matrix is represented in the generalised Gell-Mann basis of `Mat(ℂ⁴)`, composed
of Hermitian traceless matrices together with the identity matrix.
"""
ITensors.space(::SiteType"vElectron") = 16

# Shorthand notation:
function vstate(sn::StateName, ::SiteType"vElectron")
    v = ITensors.state(sn, SiteType("Electron"))
    return vec(kron(v, v'), gellmannbasis(4))
end
function vop(sn::StateName, ::SiteType"vElectron")
    return _hilbertschmidt_vec(
        op(statenamestring(sn), siteind("Electron")), gellmannbasis(4)
    )
end

# States (actual ones)
# --------------------
ITensors.state(sn::StateName"Emp", st::SiteType"vElectron") = vstate(sn, st)
ITensors.state(sn::StateName"Up", st::SiteType"vElectron") = vstate(sn, st)
ITensors.state(sn::StateName"Dn", st::SiteType"vElectron") = vstate(sn, st)
ITensors.state(sn::StateName"UpDn", st::SiteType"vElectron") = vstate(sn, st)

# States representing vectorised operators
# ----------------------------------------
ITensors.state(sn::StateName"Id", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Nup", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Ndn", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Ntot", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"NupNdn", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Aup", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Adagup", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Adn", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Adagdn", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"F", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Fup", st::SiteType"vElectron") = vop(sn, st)
ITensors.state(sn::StateName"Fdn", st::SiteType"vElectron") = vop(sn, st)

function ITensors.state(::StateName"AupF", ::SiteType"vElectron")
    return vec(
        ITensors.op(OpName("Aup"), SiteType("Electron")) *
        ITensors.op(OpName("F"), SiteType("Electron")),
        gellmannbasis(4),
    )
end
function ITensors.state(::StateName"AdagupF", ::SiteType"vElectron")
    return vec(
        ITensors.op(OpName("Adagup"), SiteType("Electron")) *
        ITensors.op(OpName("F"), SiteType("Electron")),
        gellmannbasis(4),
    )
end

function ITensors.state(::StateName"A", st::SiteType"vElectron")
    return ITensors.state(StateName("AupF"), st) + ITensors.state(StateName("Adn"), st)
end

function ITensors.state(::StateName"Adag", st::SiteType"vElectron")
    return ITensors.state(StateName("AdagupF"), st) +
           ITensors.state(StateName("Adagdn"), st)
end
