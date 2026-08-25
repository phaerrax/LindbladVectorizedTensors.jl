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

# States and operators
# ---------------------

# States derived from the Fermion site type
register_vectorized_names(
    SiteType("vFermion");
    states=("Emp", "Occ"),
    operators=("Id", "N", "F", "A", "a", "Adag", "adag", "A†", "a†"),
)

# "Up"/"Dn" are aliases of "Occ"/"Emp" for the Fermion site type.
function ITensors.state(::StateName"Up", st::SiteType"vFermion")
    return ITensors.state(StateName("Occ"), st)
end
function ITensors.state(::StateName"Dn", st::SiteType"vFermion")
    return ITensors.state(StateName("Emp"), st)
end
