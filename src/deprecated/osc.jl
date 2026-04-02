function osc_depwarnmsg()
    warnmsg = "The \"Osc\" and \"vOsc\" site types are deprecated, and will be retired in \
    a future version of the package, along with the states \"a-\", \"a+\", \"asum\" \
    and \"Asum\".  Please use the native ITensor site type \"Boson\" instead, \
    or \"vBoson\" for its vectorised counterpart."

    Base.depwarn(warnmsg, :state; force=true)
end

# Old Osc definitions
ITensors.alias(::SiteType"Osc") = SiteType"Boson"()
ITensors.alias(::SiteType"vOsc") = SiteType"vBoson"()

"""
    ITensors.space(st::SiteType"Osc";
                   dim = 2,
                   conserve_qns = false,
                   conserve_number = false,
                   qnname_number = "Number")

Create the Hilbert space for a site of type "Osc".

Optionally specify the conserved symmetries and their quantum number labels.
"""
function ITensors.space(st::SiteType"Osc"; kwargs...)
    osc_depwarnmsg()
    return ITensors.space(ITensors.alias(st); kwargs...)
end

function ITensors.val(vn::ValName, st::SiteType"Osc")
    osc_depwarnmsg()
    return ITensors.val(vn, ITensors.alias(st))
end

function ITensors.state(sn::StateName, st::SiteType"Osc", s::Index; kwargs...)
    osc_depwarnmsg()
    return ITensors.state(sn, ITensors.alias(st), s; kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"Osc", dims::Int...; kwargs...)
    osc_depwarnmsg()
    return ITensors.op(on, ITensors.alias(st), dims...; kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"Osc", s1::Index, s_tail::Index...; kwargs...)
    osc_depwarnmsg()
    return ITensors.op(on, ITensors.alias(st), s1, s_tail...; kwargs...)
end

ITensors.alias(::OpName"a-") = OpName"A"()
ITensors.alias(::OpName"a+") = OpName"Adag"()
ITensors.alias(::OpName"asum") = OpName"Asum"()

function ITensors.op(::OpName"Asum", st::SiteType"Osc", d::Int)
    osc_depwarnmsg()
    return ITensors.op(OpName("A"), st, d) + ITensors.op(OpName("Adag"), st, d)
end
function ITensors.op(on::OpName"a+", st::SiteType"Osc", d::Int)
    osc_depwarnmsg()
    return ITensors.op(ITensors.alias(on), st, d)
end
function ITensors.op(on::OpName"a-", st::SiteType"Osc", d::Int)
    osc_depwarnmsg()
    return ITensors.op(ITensors.alias(on), st, d)
end
function ITensors.op(on::OpName"asum", st::SiteType"Osc", d::Int)
    osc_depwarnmsg()
    return ITensors.op(ITensors.alias(on), st, d)
end

# Aliases from vOsc to vBoson
"""
    ITensors.space(st::SiteType"vOsc"; dim = 2)

Create the Hilbert space for a site of type "vOsc", i.e. a vectorised
harmonic oscillator of dimension `dim` (this means that the space has dimension
`dim^2`), where the vectorisation is performed wrt the (Hermitian) Gell-Mann
basis of `Mat(ℂᵈⁱᵐ)`.
"""
function ITensors.space(st::SiteType"vOsc"; kwargs...)
    osc_depwarnmsg()
    return ITensors.space(ITensors.alias(st); kwargs...)
end

function ITensors.val(vn::ValName, st::SiteType"vOsc")
    osc_depwarnmsg()
    return ITensors.val(vn, ITensors.alias(st))
end

function ITensors.state(sn::StateName, st::SiteType"vOsc", s::Index; kwargs...)
    osc_depwarnmsg()
    return ITensors.state(sn, ITensors.alias(st), s; kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"vOsc", dims::Int...; kwargs...)
    osc_depwarnmsg()
    return ITensors.op(on, ITensors.alias(st), dims...; kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"vOsc", s1::Index, s_tail::Index...; kwargs...)
    osc_depwarnmsg()
    return ITensors.op(on, ITensors.alias(st), s1, s_tail...; kwargs...)
end

export mixedlindbladplus, mixedlindbladminus

# vOsc pseudomode operators
"""
    mixedlindbladplus(n1::Int, n2::Int)

Return on OpSum expression which represents a mixed dissipation operator, appearing in the
equation for two pseudomodes, on sites at positions `n1` and `n2` in the system.
"""
function mixedlindbladplus(n1::Int, n2::Int)
    Base.depwarn(
        "This function is deprecated and will be removed in a future version.",
        :mixedlindbladplus;
        force=true,
    )
    x = OpSum()
    x += "A⋅", n1, "⋅Adag", n2
    x += "⋅Adag", n1, "A⋅", n2
    x += -0.5, "Adag⋅", n1, "A⋅", n2
    x += -0.5, "A⋅", n1, "Adag⋅", n2
    x += -0.5, "⋅Adag", n1, "⋅A", n2
    x += -0.5, "⋅A", n1, "⋅Adag", n2
    return x
end

"""
    mixedlindbladminus(n1::Int, n2::Int)

Return on OpSum expression which represents a mixed dissipation operator, appearing in the
equation for two pseudomodes, on sites at positions `n1` and `n2` in the system.
"""
function mixedlindbladminus(n1::Int, n2::Int)
    Base.depwarn(
        "This function is deprecated and will be removed in a future version.",
        :mixedlindbladminus;
        force=true,
    )
    x = OpSum()
    x += "Adag⋅", n1, "⋅A", n2
    x += "Adag⋅", n2, "⋅A", n1
    x += -0.5, "A⋅", n1, "Adag⋅", n2
    x += -0.5, "Adag⋅", n1, "A⋅", n2
    x += -0.5, "⋅A", n1, "⋅Adag", n2
    x += -0.5, "⋅Adag", n1, "⋅A", n2
    return x
end
