"""
    ITensors.space(st::SiteType"vQubit"; dim = 2)

Create the Hilbert space for a site of type "vQubit", i.e. a vectorised
qbit, where the vectorisation is performed wrt the Pauli transfer matrix
basis of `Mat(ℂ²)`, composed of the identity matrix and the three Pauli matrices.
"""
ITensors.space(::SiteType"vQubit") = 4

# States derived from the Qubit site type
register_vectorized_names(
    SiteType("vQubit"); states=("0", "1"), operators=("Id", "X", "Y", "Z", "H")
)
