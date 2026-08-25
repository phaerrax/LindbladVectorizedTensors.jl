export dissipator_loss, dissipator_gain, dissipator

"""
    ITensors.space(st::SiteType"vBoson"; dim = 2)

Create the Hilbert space for a site of type "vBoson", i.e. a vectorised
harmonic oscillator of dimension `dim` (this means that the space has dimension
`dim^2`), where the vectorisation is performed wrt the (Hermitian) Gell-Mann
basis of `Mat(ℂᵈⁱᵐ)`.
"""
function ITensors.space(::SiteType"vBoson"; dim=2)
    return dim^2
end

ITensors.op(::OpName"F", st::SiteType"vBoson", s::Index) = ITensors.op("Id", s)
# This function seems to work correctly... contrary to what I write in src/vectorize.jl,
# maybe it is possible to define operators the normal way, and not just in the
# pre/post-multiplication form...

# States
# ------

function ITensors.state(sn::StateName, st::SiteType"vBoson", s::Index; kwargs...)
    d = isqrt(ITensors.dim(s))
    stvec = state(sn, st, d; kwargs...)
    return itensor(stvec, s)
end

function ITensors.state(::StateName{N}, vst::SiteType"vBoson", d::Int) where {N}
    # Eigenstates êₙ ⊗ êₙ of the number operator, wrt the Hermitian basis.
    n = parse(Int, String(N))
    v = zeros(d)
    v[n + 1] = 1.0
    st = nonvec_stype(vst)
    return _hilbertschmidt_vec(kron(v, v'), vectorizationbasis(st, 1; dim=d))
end

# States derived from the Boson site type
register_vectorized_names(
    SiteType("vBoson"); states=(), operators=("Adag", "A", "N", "Id", "X", "Y"), dim=true
)

# Thermal states
function ITensors.state(
    ::StateName"ThermEq", vst::SiteType"vBoson", d::Int; frequency::Real, temperature::Real
)
    return if temperature == 0
        ITensors.state(StateName("0"), vst, d)
    else
        st = nonvec_stype(vst)
        numop = ITensors.op(OpName("N"), st, d)
        # We don't need to define our own matrix for the number operator when
        # we can call this one instead.
        ρ_eq = exp(-frequency / temperature * numop)
        ρ_eq /= tr(ρ_eq)
        _hilbertschmidt_vec(ρ_eq, vectorizationbasis(st, 1; dim=d))
    end
end

# Product of X = (a+a†)/√2 and of the thermal equilibrium state Z⁻¹vec(exp(-βH)).
# It is used in the computation of the correlation function of the bath.
function ITensors.state(
    ::StateName"X⋅Therm", vst::SiteType"vBoson", d::Int; frequency::Real, temperature::Real
)
    st = nonvec_stype(vst)
    xop = ITensors.op(OpName("X"), st, d)
    if temperature == 0
        ρ_eq = zeros(Float64, d, d)
        ρ_eq[1, 1] = 1.0
    else
        numop = ITensors.op(OpName("N"), st, d)
        ρ_eq = exp(-frequency / temperature * numop)
        ρ_eq /= tr(ρ_eq)
    end
    return _hilbertschmidt_vec(xop * ρ_eq, vectorizationbasis(st, 1; dim=d))
end

# GKSL equation terms
# -------------------
# For reference: given an Index `i`,
#   op(i, "B * A") == replaceprime(op(i', "B") * op(i, "A"), 2, 1)
# so the composition is performed from right to left.
# Look out for the correct order of operations:
#   vec(ABx) == op(A⋅) * vec(Bx) == op(A⋅) * (op(B⋅) * vec(x))
#   vec(xAB) == op(⋅B) * vec(xA) == op(⋅B) * (op(⋅A) * vec(x))

# Separate absorption and dissipation terms in GKSL equation
"""
    dissipator_gain(n::Int)

Return an OpSum object representing a dissipator operator (as in a GKSL equation) on site
`n`, associated to the jump operator `Adag` (the creation operator).
The OpNames `A` and `Adag` must be defined for the site type in use.
"""
function dissipator_gain(n::Int)
    dsp = OpSum()
    # a† ρ a - ½ a a† ρ - ½ ρ a a†
    dsp += "Adag⋅ * ⋅A", n
    dsp += -0.5, "A⋅ * Adag⋅", n
    dsp += -0.5, "⋅Adag * ⋅A", n
    return dsp
end

"""
    dissipator_loss(n::Int)

Return an OpSum object representing a dissipator operator (as in a GKSL equation) on site
`n`, associated to the jump operator `A` (the annihilation operator).
The OpNames `A` and `Adag` must be defined for the site type in use.
"""
function dissipator_loss(n::Int)
    dsp = OpSum()
    # a ρ a† - ½ a† a ρ - ½ ρ a† a
    dsp += "A⋅ * ⋅Adag", n
    dsp += -0.5, "N⋅", n
    dsp += -0.5, "⋅N", n
    return dsp
end

"""
    dissipator(n::Int, frequency::Real, temperature::Real)

Return an OpSum object which represents the dissipator in the GKSL equation for given
`frequency` and `temperature`, i.e. with dissipation coefficient `1/(exp(freq/temp) - 1)`.
"""
function dissipator(n::Int, frequency::Real, temperature::Real)
    # Use `expm1(x)` instead of `e^x-1` for better precision.
    avgn = 1 / expm1(frequency / temperature)
    # Julia can do math with 1/0 = Inf and 1/Inf = 0, so if temperature = 0 then avgn will
    # be 0, there's no division-by-zero error here.
    return (avgn + 1) * dissipator_loss(n) + avgn * dissipator_gain(n)
end
