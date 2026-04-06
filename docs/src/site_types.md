# Site types included in this package

This package introduces SiteTypes that allow working with MPS representing
density matrices instead of pure states.
These new SiteTypes are named by prefixing a `v` to the ITensor SiteType names.
For example, if `siteinds("Boson", N)` defines a collection of site indices to
create an MPS describing a pure state of a \\(N\\)-particle bosonic ensemble,
with `siteinds("vBoson", N)` we can create an MPS describing a _mixed_ state of
such an ensemble.
This library defines the following new SiteTypes:

- "vS=1/2"
- "vQubit"
- "vBoson"
- "vFermion"
- "vElectron"

Contrary to the native ITensor types, these ones do not support at the moment
the conservation of quantum numbers.

Internally, mixed states are represented as vectors by extracting their
coordinates with respect to the Gell-Mann basis of the appropriate dimension.
With mixed states, expectation values are given by the formula \\(\avg{A}\sb\rho
= \tr(A \rho)\\), which is also the inner product between \\(\adj{A}\\) and
\\(\rho\\) in the Hilbert space of linear operators (with respect to which the
Gell-Mann basis is an orthonormal basis). In order to be able to compute these
inner products, this package extends the meaning of “state” to include all
linear operators, so that \\(\avg{A}\sb\rho\\) is computed as the inner
product of the two “states” \\(\adj{A}\\) and \\(\rho\\).
Operators acting on states, sometimes called “superoperators”, are described in
the [Operations](@ref) section.

## "vS=1/2" SiteType

Site indices with the "vS=1/2" site type represent the density matrix of a
\\(S=1/2\\) spin.

Making a single "vS=1/2" site or collection of \\(N\\) "vS=1/2" sites:

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> s = siteind("vS=1/2");

julia> sites = siteinds("vS=1/2", N);
```

### "vS=1/2" (actual) states

The available state names for "vS=1/2" sites are:

- `"Up"` (aliases: `"↑"` and `"Z+"`) spin in the \\(\proj{\spinup}\\) state
- `"Dn"` (aliases: `"↓"` and `"Z-"`) spin in the \\(\proj{\spindown}\\) state
- `"X+"` spin in the \\(\proj{+}\\) state (projector on the +1 eigenstate of
  \\(\spin{x}\\))
- `"X-"` spin in the \\(\proj{-}\\) state (projector on the -1 eigenstate of
  \\(\spin{x}\\))
- `"Y+"` spin in the \\(\proj{\iu}\\) state (projector on the +1 eigenstate of
  \\(\spin{y}\\))
- `"Y-"` spin in the \\(\proj{-\iu}\\) state (projector on the -1 eigenstate of
  \\(\spin{y}\\))

### "vS=1/2" vectorised operators

Vectorised operators for "vS=1/2" sites can be created using the `state`
function as well, using as state name the `OpName` they would have for the
"S=1/2" site type.

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> Sx = state("Sx", s);

julia> N3 = state("N", sites[3]);
```

Available operators:

- `"Id"` Identity operator \\(I\sb2\\)
- `"σx"` (alias: `"X"`) Pauli X matrix \\(\sigma\sb{x}\\)
- `"σy"` (alias: `"Y"`) Pauli Y matrix \\(\sigma\sb{y}\\)
- `"σz"` (alias: `"Z"`) Pauli Z matrix \\(\sigma\sb{z}\\)
- `"Sx"` Spin X operator \\(\spin{x} = \frac{1}{2} \sigma\sb{x}\\)
- `"Sy"` Spin Y operator \\(\spin{y} = \frac{1}{2} \sigma\sb{y}\\)
- `"Sz"` Spin Z operator \\(\spin{z} = \frac{1}{2} \sigma\sb{z}\\)
- `"S+"` Spin raising operator \\(\spin{+} = \spin{x} + \iu\spin{y}\\)
- `"S-"` Spin lowering operator \\(\spin{-} = \spin{x} - \iu\spin{y}\\)
- `"N"` Number operator \\(N = \frac{1}{2} (I\sb2+\sigma\sb{z})\\)

## "vQubit" SiteType

A SiteType representing the mixed state of a qubit.

Making a single "vQubit" site or collection of \\(N\\) "vQubit" sites:

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> s = siteind("vQubit");

julia> sites = siteinds("vQubit", N);
```

### "vQubit" (actual) states

The available state names for "vQubit" sites are:

- `"0"` qubit in the \\(\proj{0}\\) state
- `"1"` qubit in the \\(\proj{1}\\) state

### "vQubit" vectorised operators

Vectorised operators for "vQubit" sites can be created using the `state`
function as well, using as state name the `OpName` they would have for the
"Qubit" site type.

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> X = state("X", s);

julia> H3 = state("H", sites[3]);
```

Available operators:

- `"Id"` Identity operator \\(I\sb2\\)
- `"X"` X gate
- `"Y"` Y gate
- `"Z"` Z gate
- `"H"` Hadamard gate

## "vBoson" SiteType

Site indices with the "vBoson" site type represent the density matrix of a
bosonic particle, or a harmonic oscillator (with a truncated, thus
finite-dimensional Hilbert space).

The keyword argument `dim` (default: 2) can be provided to specify the dimension
of the index, i.e. the number of available energy levels plus one.

Making a single "vBoson" site or collection of \\(N\\) "vBoson" sites:

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> s = siteind("vBoson");

julia> sites = siteinds("vBoson", N; dim=4);
```

### "vBoson" (actual) states

The available state names for "vBoson" sites are:

- `"n"` (where `n` is an integer between `0` and `dim-1`) the projector on the
  \\(n\\)-th level, \\(\proj{n}\\)
- `"ThermEq"` the thermal equilibrium state \\(\frac{1}{Z}\exp(-\frac{\omega}{T}
  N)\\); the \\(\omega\\) and \\(T\\) parameters can be given via the
  `frequency` and `temperature` keyword arguments of the `state` function,
  respectively

Examples:

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> s = siteind("vBoson"; dim=4);

julia> thermal_eq = state("ThermEq", s; temperature=1.0, frequency=5.0);

julia> fock3 = state("3", s);
```

### "vBoson" vectorised operators

Vectorised operators for "vBoson" sites can be created using the `state`
function as well, using as state name the `OpName` they would have for the
"Boson" site type.

- `"A"` annihilation operator
- `"Adag"` creation operator (beware that in the truncated space,
  \\(\adj{a}\ket{d-1}=0\\))
- `"X"` “real” quadrature \\(X=\frac{1}{\sqrt{2}}(\adj{a}+a)\\)
- `"Y"` “imaginary” quadrature \\(Y=\frac{\iu}{\sqrt{2}}(\adj{a}-a)\\)
- `"N"` number operator
- `"Id"` identity operator

## "vFermion" SiteType

Site indices with the "vFermion" SiteType represent spinless fermion sites with
the states \\(\proj0\\), \\(\proj1\\), corresponding to the empty and the
occupied state, respectively.

### "vFermion" (actual) states

Available state names for "vFermion" sites:

- `"Emp"` (alias: `"Up"`) unoccupied fermion site
- `"Occ"` (alias: `"Dn"`) occupied fermion site

### "vFermion" vectorised operators

Vectorised operators associated with "vFermion" sites:

- `"A"` (alias: `"a"`) annihilation operator
- `"Adag"` (aliases: `"adag"`, `"A†"`, `"a†"`) creation operator
- `"N"` number operator
- `"F"` Jordan-Wigner string operator
- `"Id"` identity operator

!!! warning "Anti-commutativity not implemented"
    The "vFermion" creation and annihilation operators do not account for
    fermion anticommutation (like the `"A"` and `"Adag"` Fermion operators), in
    that the Jordan-Wigner "F" operator is not automatically inserted by the
    OpSum function when building a many-body Hamiltonian operator.
    You will have to manually insert `"F"` where necessary to ensure the proper
    definition of anti-commuting operators.

## "vElectron" SiteType

The basis states of site indices with the "vElectron" SiteType correspond to
\\(\proj0\\), \\(\proj{\spinup}\\), \\(\proj{\spindown}\\) and
\\(\proj{\spinup\spindown}\\).

Making a single "vElectron" site or collection of \\(N\\) "vElectron" sites:

```jldoctest sitetypes; setup = :(using ITensors, ITensorMPS; N = 4)
julia> s = siteind("vElectron");

julia> sites = siteinds("vElectron", N);
```

### "vElectron" (actual) states

The available state names for "vElectron" sites are:

- `"Emp"` unoccupied electron site
- `"Up"` electron site occupied with one up electron
- `"Dn"` electron site occupied with one down electron
- `"UpDn"` electron site occupied with two electrons (one up, one down)

### "vElectron" vectorised operators

Vectorised operators associated with "vElectron" sites:

- `"Aup"` up-spin annihilation operator
- `"Adagup"` up-spin creation operator
- `"Adn"` down-spin annihilation operator
- `"Adagdn"` down-spin creation operator
- `"Fup"` up-spin Jordan-Wigner string operator
- `"Fdn"` down-spin Jordan-Wigner string operator
- `"F"` Jordan-Wigner string operator
- `"Nup"` number of up-spin electrons
- `"Ndn"` number of down-spin electrons
- `"Ntot"` total number operator, \\(N\sb+ + N\sb-\\)
- `"NupNdn"` product of \\(N\sb+\\) and \\(N\sb-\\)
- `"Id"` identity operator

!!! warning "Anti-commutativity not implemented"
    The "vElectron" creation and annihilation operators do not account for
    fermion anticommutation (like the `"Aup"`, `"Adn"`, `"Adagup"` and
    `"Adagdn"` Electron operators), in that the Jordan-Wigner "F" operator is
    not automatically inserted by the OpSum function when building a many-body
    Hamiltonian operator.
    You will have to manually insert `"F"` (or `"Fup"` and `"Fdn"`) where
    necessary to ensure the proper definition of anti-commuting operators.
