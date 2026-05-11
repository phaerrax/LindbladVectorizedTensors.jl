### Canonical basis

"""
    canonicalmatrix(i, j, dim)

Return the (`i`,`j`) element of the canonical basis of ``Mat(ℂᵈⁱᵐ)``, i.e. a
`dim`×`dim` matrix whose element on the `i`-th row and `j`-th column is ``1``,
and zero elsewhere.
"""
function canonicalmatrix(i, j, dim)
    m = zeros(ComplexF64, dim, dim)
    m[i, j] = 1
    return m
end

"""
    canonicalbasis(dim)

Return a list of the matrices in the canonical basis of ``Mat(ℂᵈⁱᵐ)``. 
The list is ordered corresponding to column-based vectorisation, i.e.

    canonicalbasis(dim)[j] = canonicalmatrix((j-1)%dim + 1, (j-1)÷dim + 1, dim)

with ``j ∈ {1,…,dim²}``. With this ordering,
``vec(A)ⱼ = tr(canonicalbasis(dim)[j]' * A)``.
"""
function canonicalbasis(dim; columnmajor=true)
    return if columnmajor
        [canonicalmatrix(i, j, dim) for (i, j) in [Base.product(1:dim, 1:dim)...]]
    else
        [canonicalmatrix(j, i, dim) for (i, j) in [Base.product(1:dim, 1:dim)...]]
    end
end

# Generalised Gell-Mann matrices

"""
    gellmannmatrix(j, k, dim)

Return the `(j,k)` generalised Gell-Mann matrix of dimension `dim`, normalised
wrt the Hilbert-Schmidt inner product ``(A,B) = tr(A†B)``.
The matrices are indexed as follows:

    * if ``j > k`` the matrix is symmetric and traceless;
    * if ``j < k`` the matrix is antisymmetric;
    * if ``j = k`` the matrix is diagonal.

In particular, ``j = k = dim`` gives a matrix proportional to the identity.
The two indices `j` and `k` determine the non-zero coefficients of the matrix.
The whole set of (different) Gell-Mann matrices that can be generated with this
function is a basis of ``Mat(ℂᵈⁱᵐ)``.
"""
function gellmannmatrix(j, k, dim)
    if j > dim || k > dim || j < 0 || k < 0
        throw(DomainError)
    end
    m = zeros(ComplexF64, dim, dim)
    if j > k
        m[j, k] = 1 / sqrt(2)
        m[k, j] = 1 / sqrt(2)
    elseif k > j
        m[j, k] = -im / sqrt(2)
        m[k, j] = im / sqrt(2)
    elseif j == k && j < dim
        for i in 1:j
            m[i, i] = 1
        end
        m[j + 1, j + 1] = -j
        m .*= sqrt(1 / (j * (j + 1)))
    else
        for i in 1:dim
            m[i, i] = 1 / sqrt(dim)
        end
    end
    return m
end

"""
    gellmannbasis(d, nsites=1)

Return a list containing a Hermitian basis of the tensor product of `nsites` copies of
``Mat(ℂᵈ)``, composed of (tensor products of) ``d²`` generalised Gell-Mann matrices.
"""
function gellmannbasis(dim, nsites=1)
    # Same as ptmbasis, but with a different single-site basis.
    # TODO Unify the functions in a single one that returns the multi-site basis starting
    # from the single-site one?
    single_site_basis = [
        gellmannmatrix(j, k, dim) for (j, k) in [Base.product(1:dim, 1:dim)...]
    ]

    if nsites == 1
        return single_site_basis  # and don't bother with the rest
    else
        bxn = Base.product(repeat([single_site_basis], nsites)...)
        tensorproducts = [kron(s...) for s in bxn]
        perm = reverse(ntuple(i -> i, Val{nsites}()))  # (nqbits, nqbits - 1, ..., 2, 1)
        tensorproducts_transposed = permutedims(tensorproducts, perm)

        return Base.vec(tensorproducts_transposed)
    end
end
