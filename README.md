# LindbladVectorizedTensors

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://phaerrax.github.io/LindbladVectorizedTensors.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://phaerrax.github.io/LindbladVectorizedTensors.jl/dev/)
[![Build Status](https://github.com/phaerrax/LindbladVectorizedTensors.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/phaerrax/LindbladVectorizedTensors.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package extends the [ITensors](https://github.com/ITensor/ITensors.jl)
package by defining new site types that can be used to work with systems
described by density matrices instead of pure states.

## Installation

### From a registry

This package is registered in my
[TensorNetworkSimulations](https://github.com/phaerrax/TensorNetworkSimulations)
registry. By first adding this registry, with

```julia
using Pkg
pkg"registry add https://github.com/phaerrax/TensorNetworkSimulations.git"
```

(this must be done just once per Julia installation) the package can then be
installed as a normal one:

```julia
using Pkg
Pkg.add("LindbladVectorizedTensors")
```

### From GitHub

Alternatively, straight installation from GitHub is also possible:

```julia
using Pkg
Pkg.add("https://github.com/phaerrax/LindbladVectorizedTensors.jl")
```
