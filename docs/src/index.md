# LindbladVectorizedTensors.jl

*This is the documentation for the LindbladVectorizedTensors.jl package.*

This package extends the [ITensor](https://itensor.org/) library by providing
new site types that allow working with MPS representing mixed states.
Features of this package include:

* new site types `vS=1/2`, `vBoson`, `vFermion`, `vElectron` and so on,
* functions to set up the terms in a GKSL equation easily,
* quick creation of Hamiltonian operators of common physical models.

See [Site types included in this package](@ref) for a complete list of the new
site types introduced in this package [Operations](@ref) for how to use the
vectorised types and [Examples](@ref) for some tutorials.
