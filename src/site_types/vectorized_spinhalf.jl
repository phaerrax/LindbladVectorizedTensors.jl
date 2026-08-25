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

# States derived from the S=1/2 site type
register_vectorized_names(
    SiteType("vS=1/2");
    states=("Up", "Dn", "↑", "↓", "X+", "X-", "Y+", "Y-", "Z+", "Z-"),
    operators=("Sx", "Sy", "Sz", "X", "Y", "Z", "σx", "σy", "σz", "Id", "N", "S+", "S-"),
)
