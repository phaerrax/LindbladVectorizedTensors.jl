# Operations

In this section we will see how to use the functions defined in this package to
create operators acting on mixed states.

For each of these site types, each operator from the "original", non-vectorized ITensor site
type is promoted as a pre- or post-multiplication operator.
In other words, if there exists an operator "A" for the "Fermion" type, then you will
find operators "A⋅" and "⋅A" (the dot here is a `\cdot`) available for the "vFermion"
site type.

Moreover, the `gkslcommutator` function allows you to easily define the unitary term (the
commutator) in the GKSL equation: given a list of operator names and site indices,
it returns an OpSum object representing the ``x\mapsto -i[A,x]`` map, for example
```julia-repl
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
