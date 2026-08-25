# Canonical basis

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
    gellmannbasis(d)

Return a list containing a Hermitian basis for ``Mat(ℂᵈ)`` composed of ``d²`` generalised
Gell-Mann matrices.
"""
function gellmannbasis(dim)
    return [gellmannmatrix(j, k, dim) for (j, k) in [Base.product(1:dim, 1:dim)...]]
end

# "Pauli transfer matrix" basis --- the vQubit type uses this instead of the Gell-Mann
# basis, for historical reasons.

function ptmbasis()
    st = SiteType("Qubit")
    id = Matrix(I, 2, 2)
    x = ITensors.op("X", st)
    y = ITensors.op("Y", st)
    z = ITensors.op("Z", st)
    return (1 / sqrt(2)) .* [id, x, y, z]
end

# Multi-site bases

# Transform the single-site bases above in bases for multiple sites.
function multi_site_basis(single_site_basis, nsites)
    return if nsites == 1
        single_site_basis  # nothing to do
    else
        B = Base.product(repeat([single_site_basis], nsites)...)
        # This is the Cartesian product of the single-site basis with itself `nsites` times.
        # Each element of B is a list (b_1, b_2, ..., b_n) where b_i is a basis matrix.
        # For example, with two sites, if
        #   tp(a, b) = "$a ⊗ $b"
        # and
        #   ssb = ["s$n" for n in 0:3]
        # then we obtain
        #   [tp(b...) for b in Base.product(repeat([ssb], 2)...)] =
        #    "b0 ⊗ b0"  "b0 ⊗ b1"  "b0 ⊗ b2"  "b0 ⊗ b3"
        #    "b1 ⊗ b0"  "b1 ⊗ b1"  "b1 ⊗ b2"  "b1 ⊗ b3"
        #    "b2 ⊗ b0"  "b2 ⊗ b1"  "b2 ⊗ b2"  "b2 ⊗ b3"
        #    "b3 ⊗ b0"  "b3 ⊗ b1"  "b3 ⊗ b2"  "b3 ⊗ b3"

        # By calling `Base.vec` on this matrix we stack its columns, but we want to unroll
        # by rows instead, so we transpose it first. The call to `permutedims` below
        # does this transposition "on all dimensions".
        # (Why do we need rows instead of columns here?)
        tensorproducts = [kron(b...) for b in B]
        perm = reverse(ntuple(i -> i, Val(nsites)))  # (nsites, nsites - 1, ..., 2, 1)
        tensorproducts_transposed = permutedims(tensorproducts, perm)

        # We transform each of them in the tensor product b_1 ⊗ b_2 ⊗ ... ⊗ b_n.
        Base.vec(tensorproducts_transposed)
    end
end

# All methods are defined with the `dim` keyword argument, for simplicity reasons, even if
# makes sense for the Boson type only (the others use the default correct value).
function vectorizationbasis(st::SiteType, nsites::Int; dim=nothing)
    return multi_site_basis(gellmannbasis(ITensors.space(st)), nsites)
end

# Explicit override for Bosons, which need the site dimension supplied.
function vectorizationbasis(st::SiteType"Boson", nsites::Int; dim)
    return multi_site_basis(gellmannbasis(ITensors.space(st; dim)), nsites)
end

# Explicit override for Qubits, which use the PTM basis instead.
function vectorizationbasis(st::SiteType"Qubit", nsites::Int; dim=nothing)
    return multi_site_basis(ptmbasis(), nsites)
end
