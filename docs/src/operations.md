# Operations

In this section we will see how to use the functions defined in this package to
create operators acting on mixed states.

## Left- and right-multiplication

For each vectorised site type, each operator from the “original”, non-vectorised
site type can be used to multiply a state (in the sense of this package, that is
a mixed state or an operator) on the left on the right.  In other words, if
there exists an operator `"A"` for the "T" type, then the operators `"A⋅"` and
`"⋅A"` (the dot is a `\cdot`) are automatically defined for the "vT" site
type, which perform left- and right-multiplication by `"A"`, respectively.

For example, let's play around with a vFermion site. We create a site index and
define on it the empty state \\(\proj{0}\\), then apply the creation operator on
the left.

```jldoctest operations; setup = :(using LindbladVectorizedTensors, ITensors, ITensorMPS)
julia> s = siteind("vFermion");

julia> emp = state(s, "Emp");

julia> Adagemp = apply(op("Adag⋅", s), emp);

```

We have now obtained \\(\adj{a}\proj{0}\\), and we can check for example that its
inner product with \\(\proj1\\) is zero, that is, \\(\tr(\proj1 \adj{a}\proj{0})
= \braket{0}{1} \bra1 \adj{a}\ket0 = 0\\):

```jldoctest operations
julia> occ = state(s, "Occ");

julia> scalar(occ * Adagemp)
0.0 + 0.0im

```

If we multiply \\(\adj{a}\proj{0}\\) by \\(a\\) on the right too, then we get
\\(\tr(\proj1 \adj{a}\proj{0}a) = \bra0 a \ket1 \bra1 \adj{a}\ket0 = 1\\):

```jldoctest operations
julia> AdagempA = apply(op("⋅A", s), Adagemp);

julia> scalar(occ * AdagempA)
0.9999999999999998 + 0.0im
```

Left- and right-multiplication operators can be composed with the usual ITensor
rules: for example, they can be multiplied or added together.
In fact, we could have defined `AdagempA` in one shot by applying the
`"Adag⋅ * ⋅A"` operator:

```jldoctest operations
julia> AdagempA = apply(op("Adag⋅ * ⋅A", s), emp);

```

These, together with the identity `"Id"`, are the only “basic” operators defined
on the vectorised types which can be called with the `op` function.

!!! info "Defining new operators"
    Due to technical limitations, it is not possible to define new operators by
    overloading the `op` method for the vectorised site types in this package.

## Commutators

The time evolution of mixed states is often described by the von Neumann
equation, or the GKSL equation, both of which involve the commutator of the
state and an operator.
For this reason, this package provides a useful feature with which we can create
the \\(\rho \mapsto -\iu [A,\rho]\\) “superoperator”: this is the
`gkslcommutator` function.

Given a list of operator names and integers, in an alternating fashion like
`"A1, n1, A2, n2, ..."`, the `gkslcommutator` function returns an OpSum object
that represents the \\(-\iu[A,\blank]\\) map where \\(A\\) is the product of
`A1` on site `n1`, `A2` on site `n2` and so on, in a similar syntax as OpSum
itself.  For example,

```jldoctest operations
julia> gkslcommutator("A", 1)
sum(
  0.0 - 1.0im A⋅(1,)
  0.0 + 1.0im ⋅A(1,)
)

julia> gkslcommutator("A", 1, "B", 3)
sum(
  0.0 - 1.0im A⋅(1,) B⋅(3,)
  0.0 + 1.0im ⋅A(1,) ⋅B(3,)
)

```

These OpSum objects can then be turned into MPOs with the usual syntax
`MPO(opsum, sites)`.

If we also provide an array of site indices, we can use the
`gkslcommutator_itensor` function, that directly gives us an ITensor object.

```jldoctest operations; filter = r"id=\d+" => "id=#"
julia> sites = siteinds("vS=1/2", 3);

julia> gkslcommutator_itensor(sites, "Sx", 2)
ITensor ord=2 (dim=4|id=652|"Site,n=2,vS=1/2")' (dim=4|id=652|"Site,n=2,vS=1/2")
NDTensors.Dense{ComplexF64, Vector{ComplexF64}}

```

### Example: GKSL equation

Consider the Hamiltonian operator

```math
H = \sum_{n=1}^{N} \omega_n\phantomadj \adj{a_n} a_n\phantomadj
```

for a harmonic oscillator. Let's choose 8 as the dimension of the Hilbert space,
\\(N=3\\) and all \\(\omega\sb{n}\\)s equal to one, for simplicity. The
\\(-\iu[H,\blank]\\) map in the von Neumann equation \\(\dot\rho\sb{t} =
-\iu[H,\rho\sb{t}]\\) can be defined as follows.

```jldoctest operations
julia> N = 3; s = siteinds("vBoson", N; dim=8); ω = ones(N);

julia> h = OpSum();

julia> for n in 1:N
           h += ω[n] * gkslcommutator("N", n)
       end

julia> h
sum(
  0.0 - 1.0im N⋅(1,)
  0.0 + 1.0im ⋅N(1,)
  0.0 - 1.0im N⋅(2,)
  0.0 + 1.0im ⋅N(2,)
  0.0 - 1.0im N⋅(3,)
  0.0 + 1.0im ⋅N(3,)
)

```

We could also add to the Hamiltonian an exchange interaction

```math
\sum_{n=1}^{N} \lambda_n\phantomadj (\adj{a_n} a_{n+1}\phantomadj + \adj{a_{n+1}}
a_n\phantomadj)
```

with the following commands (we set \\(\lambda\sb{n} = \frac12\\)).

```jldoctest operations
julia> λ = fill(0.5, N);

julia> for n in 1:N-1
           h += λ[n] * gkslcommutator("Adag", n, "A", n+1)
           h += λ[n] * gkslcommutator("A", n, "Adag", n+1)
       end

julia> h
sum(
  0.0 - 1.0im N⋅(1,)
  0.0 + 1.0im ⋅N(1,)
  0.0 - 1.0im N⋅(2,)
  0.0 + 1.0im ⋅N(2,)
  0.0 - 1.0im N⋅(3,)
  0.0 + 1.0im ⋅N(3,)
  0.0 - 0.5im Adag⋅(1,) A⋅(2,)
  0.0 + 0.5im ⋅Adag(1,) ⋅A(2,)
  0.0 - 0.5im A⋅(1,) Adag⋅(2,)
  0.0 + 0.5im ⋅A(1,) ⋅Adag(2,)
  0.0 - 0.5im Adag⋅(2,) A⋅(3,)
  0.0 + 0.5im ⋅Adag(2,) ⋅A(3,)
  0.0 - 0.5im A⋅(2,) Adag⋅(3,)
  0.0 + 0.5im ⋅A(2,) ⋅Adag(3,)
)

```

## The adjoint map

Another common operation on mixed states is \\(\rho \mapsto X \rho \adj{X}\\),
where \\(X\\) is some operator. For example, in the definition of the
dissipation terms in the GKSL equation (see the [relative example](@ref "GKSL
equation") we have \\(a\sb{n}\phantomadj \rho \adj{a\sb{n}}\\), which can be
simply written as `apply(op("A⋅ * ⋅Adag", s, n), ρ)`, where `s` is the list of
Index object that `ρ` is defined on.
This was straightforward because the `A` and `Adag` operators are already
defined, but not all operators have their adjoint available in the library.

The `adjointmap_itensor` [^adjointname] function provides a convenient way of
creating the tensor representing the above map for any operator that can be
called on the non-vectorised site type.  The syntax is similar to the one of the
`op` function: for example, we can write

```jldoctest operations
julia> s = siteinds("vQubit", 6);

julia> ρ = MPS(s, "0");

julia> apply(adjointmap_itensor("CNOT", s[1], s[2]), ρ);

julia> apply(adjointmap_itensor("Rx", s, 5; θ=1/2), ρ);

julia> apply(adjointmap_itensor("Rx", s, 5; θ=pi/3), ρ);

```

We did not need to define the adjoint of "Rx": the `adjointmap_itensor` already
takes care of computing it internally, without defining a new operator.
The function automatically computes the adjoint of the given operator and
creates the tensor corresponding to the vectorisation of the \\(\rho \mapsto X
\rho \adj{X}\\) map, that can then be applied to MPSs or other ITensors.

!!! warning "Only available for the `vQubit` site tipe"
    At the time of writing this, the `adjointmap_itensor` function is available
    only for `vQubit` sites.

[^adjointname]: The name of this function comes from the similarity of this map
to the [adjoint representation]
(https://en.wikipedia.org/wiki/Adjoint_representation) of a group (typically a
Lie group).
