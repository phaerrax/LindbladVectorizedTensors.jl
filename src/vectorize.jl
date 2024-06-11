"""
    vec_purestate_densitymatrix(::SiteType"Qubit", x::MPS; existing_sites=nothing)

Return the MPS representing the vectorization of `projection(x)`, i.e. of the density matrix
``|x⟩⟨x|``, in the PTM basis.

The returned MPS will have a new set of indices; it may be desirable instead to "anchor"
the result to an already existing set of indices, i.e. for comparing it with some other
states; the `existing_sites` keyword arguments allows the caller to specify the
set of vectorized site indices that the result should have.
"""
function vec_purestate_densitymatrix(::SiteType"Qubit", x::MPS; existing_sites=nothing)
    N = length(x)
    rho = outer(x', x)  # Density matrix of the pure state `x`
    sites = siteinds(x)

    # The "combiner" index merges the two matrix indices in a single one, effectively
    # transforming the matrix into a vector. Our vQubit site type works in the PTM basis
    # though... we need to implement a change of basis on each site.
    #
    # If t = [a b; c d] with index `s` then
    #   t * combiner(s', s) = [a, c, b, d]  (stacks columns)
    #   t * combiner(s, s') = [a, b, c, d]  (stacks rows)
    #
    # Let's stack by rows, and define the canonical basis as
    #   e[1] = [1 0; 0 0]
    #   e[2] = [0 1; 0 0]
    #   e[3] = [0 0; 1 0]
    #   e[4] = [0 0; 0 1]
    # so that [a b; c d] = a*e[1] + b*e[2] + c*e[3] + d*e[4].
    # The PTM basis (at least, our version of it) is
    #   f[1] = 1/sqrt(2) [1 0; 0 1]
    #   f[2] = 1/sqrt(2) [0 1; 1 0]
    #   f[3] = 1/sqrt(2) [0 -i; i 0]
    #   f[4] = 1/sqrt(2) [1 0; 0 -1]
    # in other words, f[1] = 1/√2 I_2 and f[i+1] = 1/√2 σ_i for i∈{1,2,3}.
    # and with these definitions we get
    #   e[i] = sum([U[i,j] * f[j] for j in 1:4])
    # with
    #        1  ⎛ 1  0  0  1 ⎞
    #   U = ⸺   ⎜ 0  1  i  0 ⎟
    #       √2  ⎜ 0  1 -i  0 ⎟
    #           ⎝ 1  0  0 -1 ⎠
    # Note that coordinate vectors change basis with the transpose of this matrix,
    # so we need to multiply each blocks of the MPS by U on the _left_.
    # We achieve this by defining the tensor with a specific order of the indices.
    x_vec = MPS(length(x))  # undef MPS
    sites_rho = noprime(siteinds(first, rho))
    # `rho` is an MPO, so with `siteinds(rho)` we get both each "real" site index and its
    # primed copy. With `first` we select only one of them; I don't know if we get s or s'
    # consistently, so we use `noprime` for safety.
    for i in 1:length(x_vec)
        x_vec[i] = rho[i] * combiner(sites_rho[i]', sites_rho[i])
    end
    changeofbasis_mat = 1 / sqrt(2) .* [
        1 0 0 1
        0 1 im 0
        0 1 -im 0
        1 0 0 -1
    ]
    # The combiner produces an Index of type `(dim=4|id=...|"CMB,Link")`, but we
    # need a "vQubit" Index on each site instead. We need to change that.
    sites_combined = siteinds(x_vec)  # These will be "CMB" indices
    if isnothing(existing_sites)
        # Create a new set of indices
        sites_vec = siteinds("vQubit", length(x_vec))
    else
        if length(x_vec) != length(existing_sites)
            error(
                "Provided sites for the vectorized state are not compatible with the new MPS",
            )
        end
        sites_vec = existing_sites
    end

    changeofbasis = [
        ITensor(ComplexF64, changeofbasis_mat, sites_combined[i], sites_vec[i]) for
        i in 1:length(sites_vec)
    ]
    # TODO Rename site labels so the "n=#" part is the same as the original MPS?
    for i in 1:length(x_vec)
        x_vec[i] *= changeofbasis[i]
    end

    return x_vec
end
