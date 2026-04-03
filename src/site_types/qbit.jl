# Additional Qubit gates
# ----------------------
ITensors.op(::OpName"Id", ::SiteType"Qubit") = [
    1 0
    0 1
]
ITensors.op(::OpName"id", st::SiteType"Qubit") = ITensors.op(OpName("Id"), st)

function ITensors.op(::OpName"CH", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("Id"), st)
    h = op(OpName("H"), st)
    return kron(proj0, id) + kron(proj1, h)
end

function ITensors.op(::OpName"CPhase", st::SiteType"Qubit"; ϕ::Real)
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("Id"), st)
    phase = op(OpName("Phase"), st; ϕ=ϕ)
    return kron(proj0, id) + kron(proj1, phase)
end

function ITensors.op(::OpName"CU3", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    u = ITensors.op(OpName("U"), st; θ=θ, ϕ=ϕ, λ=λ)
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("id"), st)
    return kron(proj0, id) + kron(proj1, u)
end

function ITensors.op(::OpName"CU", st::SiteType"Qubit"; θ::Real, ϕ::Real, λ::Real)
    return ITensors.op(OpName("CU3"), st; θ=θ, ϕ=ϕ, λ=λ)
end

function ITensors.op(::OpName"CU1", st::SiteType"Qubit"; λ::Real)
    return ITensors.op(OpName("CU3"), st; θ=0, ϕ=0, λ=λ)
end

function ITensors.op(::OpName"CCCCNOT", st::SiteType"Qubit")
    proj0 = op(OpName("Proj0"), st)
    proj1 = op(OpName("Proj1"), st)
    id = op(OpName("Id"), st)
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
