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

# States and operators
# ---------------------

# States derived from the Boson site type
register_vectorized_names(
    SiteType("vElectron");
    states=("Emp", "Up", "Dn", "UpDn"),
    operators=(
        "Id",
        "Nup",
        "Ndn",
        "Ntot",
        "NupNdn",
        "Aup",
        "Adagup",
        "Adn",
        "Adagdn",
        "F",
        "Fup",
        "Fdn",
    ),
)

function ITensors.state(::StateName"AupF", ::SiteType"vElectron")
    return _hilbertschmidt_vec(
        ITensors.op(OpName("Aup"), SiteType("Electron")) *
        ITensors.op(OpName("F"), SiteType("Electron")),
        gellmannbasis(4),
    )
end
function ITensors.state(::StateName"AdagupF", ::SiteType"vElectron")
    return _hilbertschmidt_vec(
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
