export vec_projector

function _change_of_basis_to_gellmann(new_s::Index, old_s::Index)
    @assert ITensors.dim(new_s) == ITensors.dim(old_s)
    d = isqrt(ITensors.dim(new_s))

    # Let {e[n]}ₙ be the canonical (row-major) basis of d×d matrices, i.e. for d=2
    #
    #   e[1] = [1 0; 0 0]
    #   e[2] = [0 1; 0 0]
    #   e[3] = [0 0; 1 0]
    #   e[4] = [0 0; 0 1]
    #
    # so that [a b; c d] = a*e[1] + b*e[2] + c*e[3] + d*e[4].

    eb = canonicalbasis(d; columnmajor=false)

    # The Gell-Mann matrix basis is
    #
    #   f[1] = 1/sqrt(2) [0 1; 1 0]
    #   f[2] = 1/sqrt(2) [0 -i; i 0]
    #   f[3] = 1/sqrt(2) [1 0; 0 -1]
    #   f[4] = 1/sqrt(2) [1 0; 0 1],

    gb = gellmannbasis(d)

    # To go from the canonical to the Gell-Mann basis we use (again with d=2)
    #
    #   ⎛ f[1] ⎞    1  ⎛ 0  1  1  0 ⎞ ⎛ e[1] ⎞
    #   ⎜ f[2] ⎟ = ⸺   ⎜ 0 -i  i  0 ⎟ ⎜ e[2] ⎟
    #   ⎜ f[3] ⎟   √2  ⎜ 1  0  0 -1 ⎟ ⎜ e[3] ⎟.
    #   ⎝ f[4] ⎠       ⎝ 1  0  0  1 ⎠ ⎝ e[4] ⎠
    #
    # Note that the n-th row of the matrix is made up of the stacked rows of f[n].
    # We can easily generalise this equation, then, to
    #
    #   ⎛ f[1]  ⎞    ⎛ ⟨e[1],f[1]⟩  ⟨e[2],f[1]⟩  ⋯ ⟨e[d²],f[1]⟩  ⎞ ⎛ e[1]  ⎞
    #   ⎜ f[2]  ⎟ =  ⎜ ⟨e[1],f[2]⟩  ⟨e[2],f[2]⟩  ⋯ ⟨e[d²],f[2]⟩  ⎟ ⎜ e[2]  ⎟
    #   ⎜  ⋮    ⎟    ⎜      ⋮            ⋮       ⋱       ⋮       ⎟ ⎜  ⋮    ⎟.
    #   ⎝ f[d²] ⎠    ⎝ ⟨e[1],f[d²]⟩ ⟨e[2],f[d²]⟩ ⋯ ⟨e[d²],f[d²]⟩ ⎠ ⎝ e[d²] ⎠

    # Now, if ve and vf are the coordinate vectors of v in the {e[n]}ₙ and {f[n]}ₙ bases,
    # respectively, then
    #
    #   ve[k] = ∑ₗ U[l,k] vf[l],   i.e.,   ve = Uᵀ vf.
    #
    # We want vf, but the combiner gives us ve, thus we need to compute
    #
    #   vf = (Uᵀ)* ve = Ū ve.

    u = ITensor(ComplexF64, new_s, old_s)
    for i in 1:(d ^ 2), j in 1:(d ^ 2)
        u[new_s => i, old_s => j] = conj(tr(eb[j]' * gb[i]))
        # When we construct an `op` like this through a matrix, the primed index s' runs
        # over the rows while s runs over the columns, i.e.
        #
        #   m = op(s, [a b; c d])
        #   m[s' => 1, s => 1] = a
        #   m[s' => 1, s => 2] = b
        #   m[s' => 2, s => 1] = c
        #   m[s' => 2, s => 2] = d
    end

    return u
end

function _change_of_basis_to_ptm(new_s::Index, old_s::Index)
    # Like _change_of_basis_matrix_canonical, but with the PTM basis instead.
    d = 2
    eb = canonicalbasis(d)
    gb = ptmbasis(1)

    u = ITensor(ComplexF64, new_s, old_s)
    for i in 1:(d ^ 2), j in 1:(d ^ 2)
        u[new_s => i, old_s => j] = conj(tr(eb[j]' * gb[i]))
    end

    return u
end

_change_of_basis_matrix(::SiteType, ::Index, ::Index) = nothing
function _change_of_basis_matrix(::SiteType"Boson", new_s::Index, old_s::Index)
    _change_of_basis_to_gellmann(new_s, old_s)
end
function _change_of_basis_matrix(::SiteType"Fermion", new_s::Index, old_s::Index)
    _change_of_basis_to_gellmann(new_s, old_s)
end
function _change_of_basis_matrix(::SiteType"S=1/2", new_s::Index, old_s::Index)
    _change_of_basis_to_gellmann(new_s, old_s)
end
function _change_of_basis_matrix(::SiteType"Electron", new_s::Index, old_s::Index)
    _change_of_basis_to_gellmann(new_s, old_s)
end
function _change_of_basis_matrix(::SiteType"FDot3", new_s::Index, old_s::Index)
    _change_of_basis_to_gellmann(new_s, old_s)
end
function _change_of_basis_matrix(::SiteType"Qubit", new_s::Index, old_s::Index)
    _change_of_basis_to_ptm(new_s, old_s)
end

const _implemented_vtypes = ["Boson", "Fermion", "S=1/2", "Electron", "FDot3", "Qubit"]

"""
    vec_projector(x::MPS; existing_sites=nothing, kwargs...)

Return the MPS representing the vectorisation of `projector(x)`, i.e. of the density matrix
``|x⟩⟨x|``.  This method is available for the following site types:
$(join(_implemented_vtypes, ", ", " and ")).  The vectorisation is performed using the
Gell-Mann basis for all types except Qubit, for which the Pauli transfer-matrix basis is
used instead. In any case, the basis is always the same one used in the definition of the
vectorised states and operators.

The returned MPS will have a brand new set of indices, of the corresponding vectorised site
type (e.g. if `x` has Qubit indices then the result will have vQubit indices). Set
`existing_sites` to an existing set of site indices in order to use them in the resulting
MPS.

Use additional keyword arguments to control the level of truncation, which are the same as
those accepted by `contract(::MPO, ::MPO; kw...)`.

# Keyword arguments

* `existing_sites`: replace the indices of the output MPS to a desired existing set.
* `normalize::Bool=true`: whether or not to normalise the input MPS before forming the
  projector. If `normalize==false` and the input MPS is not already normalised, this
  function will not output a proper projector, and simply outputs the vectorisation of
  ``|x⟩⟨x|``, i.e. the projector scaled by `norm(x)^2`.
* truncation keyword arguments accepted by `contract(::MPO, ::MPO; kw...)`.

"""
function vec_projector(x::MPS; existing_sites=nothing, projector_kwargs...)
    N = length(x)
    projx = projector(x; projector_kwargs...)  # Density matrix of the pure state `x`
    s = siteinds(x)

    # `projx` is an MPO, so with `siteinds(projx)` we get both each "real" site index and
    # its primed copy. With `first` we select only one of them; I don't know if we get s or
    # s', so let's use `noprime` just to be sure.
    s_vec = noprime(siteinds(first, projx))

    # We want to create an MPS of the vectorised site type, which on site n contains
    # the vector [⟨f[1], x[n]⟩, ⟨f[2], x[n]⟩, ..., ⟨f[d²], x[n]⟩] where {f[i]}ᵢ is the
    # Gell-Mann basis of the space of d×d matrices.
    #
    # Each x[n] tensor has two physical (site) indices, s[n] and s[n'], and the first thing
    # we can do is to combine them into a single one. We use ITensor's "combiner".
    # The "combiner" index merges the two matrix indices into a single one, effectively
    # transforming the matrix into a vector.
    #
    # For example, assume we work in dimension 2, and let t = op(sᵢ, [a b; c d]), where sᵢ
    # is a site index.  Then
    #   t * combiner(sᵢ', sᵢ) = [a, c, b, d]  (stacks columns)
    #   t * combiner(sᵢ, sᵢ') = [a, b, c, d]  (stacks rows)
    #
    # We choose the second form of the combiner to combine the MPS indices.
    x_vec = MPS(N)  # empty (undef) MPS
    for i in 1:length(x_vec)
        x_vec[i] = projx[i] * combiner(s_vec[i], s_vec[i]')
    end

    # The combiner produces an Index of type `(dim=...|id=...|"CMB,Link")`, but we
    # need a "vSomething" index on each site instead.
    s_cmb = siteinds(x_vec)  # These will be "CMB" indices

    # Determine the SiteType of each site of the input MPS.
    # We match the found SiteTypes against the types supported by this package: there must
    # be only one result.
    # - More than one result: I don't see how it could happen if one uses ITensors
    #   functions without doing anything weird...
    # - Zero results: the original site type is not supported by this package.
    sitetypes = Vector{String}(undef, N)
    for n in 1:N
        sitetypenames_thissite = sitetypestring.(ITensors.SiteTypes._sitetypes(s[n]))
        filter!(in(_implemented_vtypes), sitetypenames_thissite)
        sitetypes[n] = only(sitetypenames_thissite)
    end

    # Create a new set of indices, or use the one provided in the keyword argument.
    s_vec = if isnothing(existing_sites)
        [
            if dim(s[n]) > 2
                siteind("v" * sitetypes[n]; dim=dim(s[n]))
            else
                siteind("v" * sitetypes[n])
            end for n in 1:N
        ]
    else
        if N != length(existing_sites)
            throw(ArgumentError("provided sites for the vectorized state are not \
                                compatible with the new MPS"))
        end
        existing_sites
    end

    u = Vector{ITensor}(undef, N)
    for n in 1:N
        cbm = _change_of_basis_matrix(SiteType(sitetypes[n]), s_vec[n], s_cmb[n])
        if !isnothing(cbm)
            u[n] = cbm
        else
            throw(
                ArgumentError(
                    "Change-of-basis matrix not found for Index tags: $(tags.(s_vec[n]))."
                ),
            )
        end
    end

    # TODO Rename site labels so the "n=#" part is the same as the original MPS?
    for n in 1:N
        x_vec[n] *= u[n]
    end

    return x_vec
end
