# Additional Qubit gates
# ----------------------

# How to combine indices of multi-site operators:
#
#   julia> proj0 = op(OpName("Proj0"), st);
#
#   julia> proj1 = op(OpName("Proj1"), st);
#
#   julia> id = matrix(op("Id", siteind(st)));
#
#   julia> r = matrix(op("Ry", s, 3; θ=pi/4))
#   2×2 Matrix{Float64}:
#    0.92388   -0.382683
#    0.382683   0.92388
#
#   julia> C = combiner(inds(cp; plev=0))
#   ITensor ord=3 (dim=4|id=48|"CMB,Link") (dim=2|id=531|"Qubit,Site,n=3")
#       (dim=2|id=895|"Qubit,Site,n=2")
#   NDTensors.Combiner
#
#   julia> matrix(C' * cr * C)
#   4×4 Matrix{Float64}:
#    1.0  0.0  0.0        0.0
#    0.0  1.0  0.0        0.0
#    0.0  0.0  0.92388   -0.382683
#    0.0  0.0  0.382683   0.92388
#
#   julia> kron(proj0, id) + kron(proj1, r)
#   4×4 Matrix{Float64}:
#    1.0  0.0  0.0        0.0
#    0.0  1.0  0.0        0.0
#    0.0  0.0  0.92388   -0.382683
#    0.0  0.0  0.382683   0.92388
#
# but matrix(C * cr * C') != kron(proj0, id) + kron(proj1, r)

function ITensors.op(::OpName"CH", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = matrix(op("Id", siteind(st)))
    h = op(OpName("H"), st)
    return kron(proj0, id) + kron(proj1, h)
end

ITensors.op(::OpName"CPhase", st::SiteType"Qubit"; kwargs...) = op("CPHASE", st; kwargs...)

function ITensors.op(::OpName"CCCCNOT", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = matrix(op("Id", siteind(st)))
    not = op(OpName("X"), st)
    return kron(proj0, proj0, proj0, proj0, id) +
           kron(proj0, proj0, proj0, proj1, id) +
           kron(proj0, proj0, proj1, proj0, id) +
           kron(proj0, proj0, proj1, proj1, id) +
           kron(proj0, proj1, proj0, proj0, id) +
           kron(proj0, proj1, proj0, proj1, id) +
           kron(proj0, proj1, proj1, proj0, id) +
           kron(proj0, proj1, proj1, proj1, id) +
           kron(proj1, proj0, proj0, proj0, id) +
           kron(proj1, proj0, proj0, proj1, id) +
           kron(proj1, proj0, proj1, proj0, id) +
           kron(proj1, proj0, proj1, proj1, id) +
           kron(proj1, proj1, proj0, proj0, id) +
           kron(proj1, proj1, proj0, proj1, id) +
           kron(proj1, proj1, proj1, proj0, id) +
           kron(proj1, proj1, proj1, proj1, not)
end
