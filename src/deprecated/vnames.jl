function vnames_depwarnmsg()
    warnmsg = "For vectorised site types, state names beginning with `v` or `vec` are \
    deprecated and will be removed in a future version of this package. States \
    corresponding to operators of the corresponding non-vectorised site type will use the \
    same name, without `v` or `vec` added as a prefix. For example, the vectorised version \
    of the operator \"Sx\" for the \"S=1/2\" type will use the StateName \"Sx\" for the \
    \"vS=1/2\" type, instead of the current \"vSx\"."
    Base.depwarn(warnmsg, :state; force=true)
end

# vElectron
function ITensors.state(::StateName"vId", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
function ITensors.state(::StateName"vNup", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Nup"), st)
end
function ITensors.state(::StateName"vNdn", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Ndn"), st)
end
function ITensors.state(::StateName"vNtot", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Ntot"), st)
end
function ITensors.state(::StateName"vNupNdn", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("NupNdn"), st)
end
function ITensors.state(::StateName"vAup", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Aup"), st)
end
function ITensors.state(::StateName"vAdagup", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adagup"), st)
end
function ITensors.state(::StateName"vAdn", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adn"), st)
end
function ITensors.state(::StateName"vAdagdn", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adagdn"), st)
end
function ITensors.state(::StateName"vF", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("F"), st)
end
function ITensors.state(::StateName"vAupF", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("AupF"), st)
end
function ITensors.state(::StateName"vAdagupF", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("AdagupF"), st)
end
function ITensors.state(::StateName"vA", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("A"), st)
end
function ITensors.state(::StateName"vAdag", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adag"), st)
end
function ITensors.state(::StateName"vecId", st::SiteType"vElectron")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end

# vFdot3
function ITensors.state(::StateName"vId", st::SiteType"vFDot3")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
function ITensors.state(::StateName"vecId", st::SiteType"vFDot3")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
function ITensors.state(::StateName"vn1", st::SiteType"vFDot3")
    vnames_depwarnmsg()
    return ITensors.state(StateName("n1"), st)
end
function ITensors.state(::StateName"vn2", st::SiteType"vFDot3")
    vnames_depwarnmsg()
    return ITensors.state(StateName("n2"), st)
end
function ITensors.state(::StateName"vn3", st::SiteType"vFDot3")
    vnames_depwarnmsg()
    return ITensors.state(StateName("n3"), st)
end
function ITensors.state(::StateName"vntot", st::SiteType"vFDot3")
    vnames_depwarnmsg()
    return ITensors.state(StateName("ntot"), st)
end

# vFermion
function ITensors.state(::StateName"vId", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
function ITensors.state(::StateName"vN", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("N"), st)
end
function ITensors.state(::StateName"vF", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("F"), st)
end
function ITensors.state(::StateName"vA", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("A"), st)
end
function ITensors.state(::StateName"va", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("a"), st)
end
function ITensors.state(::StateName"vAdag", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adag"), st)
end
function ITensors.state(::StateName"vadag", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("adag"), st)
end
function ITensors.state(::StateName"vA†", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("A†"), st)
end
function ITensors.state(::StateName"va†", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("a†"), st)
end
function ITensors.state(::StateName"vecId", st::SiteType"vFermion")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end

# vOsc
function ITensors.state(::StateName"vAdag", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adag"), st, d)
end
function ITensors.state(::StateName"vA", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("A"), st, d)
end
function ITensors.state(::StateName"vN", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("N"), st, d)
end
function ITensors.state(::StateName"vId", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st, d)
end
function ITensors.state(::StateName"vX", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("X"), st, d)
end
function ITensors.state(::StateName"vY", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("Y"), st, d)
end
function ITensors.state(::StateName"veca+", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adag"), st, d)
end
function ITensors.state(::StateName"veca-", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("A"), st, d)
end
function ITensors.state(::StateName"vecplus", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("Adag"), st, d)
end
function ITensors.state(::StateName"vecminus", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("A"), st, d)
end
function ITensors.state(::StateName"vecN", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("N"), st, d)
end
function ITensors.state(::StateName"vecId", st::SiteType"vOsc", d::Int)
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st, d)
end

ITensors.alias(::SiteType"HvOsc") = SiteType"vOsc"()
ITensors.alias(::SiteType"vecOsc") = SiteType"vOsc"()

function ITensors.space(st::SiteType"HvOsc"; kwargs...)
    vecHv_depwarnmsg("Osc")
    return ITensors.space(ITensors.alias(st); kwargs...)
end
function ITensors.space(st::SiteType"vecOsc"; kwargs...)
    vecHv_depwarnmsg("Osc")
    return ITensors.space(ITensors.alias(st); kwargs...)
end

function ITensors.val(vn::ValName, st::SiteType"HvOsc")
    vecHv_depwarnmsg("Osc")
    return ITensors.val(vn, ITensors.alias(st))
end
function ITensors.val(vn::ValName, st::SiteType"vecOsc")
    vecHv_depwarnmsg("Osc")
    return ITensors.val(vn, ITensors.alias(st))
end

function ITensors.state(sn::StateName, st::SiteType"vecOsc", s::Index; kwargs...)
    vecHv_depwarnmsg("Osc")
    return ITensors.state(sn, ITensors.alias(st), s; kwargs...)
end
function ITensors.state(sn::StateName, st::SiteType"HvOsc", s::Index; kwargs...)
    vecHv_depwarnmsg("Osc")
    return ITensors.state(sn, ITensors.alias(st), s; kwargs...)
end

function ITensors.op(
    on::OpName, st::SiteType"vecOsc", s1::Index, s_tail::Index...; kwargs...
)
    vecHv_depwarnmsg("Osc")
    return ITensors.op(on, ITensors.alias(st), s1, s_tail...; kwargs...)
end
function ITensors.op(
    on::OpName, st::SiteType"HvOsc", s1::Index, s_tail::Index...; kwargs...
)
    vecHv_depwarnmsg("Osc")
    return ITensors.op(on, ITensors.alias(st), s1, s_tail...; kwargs...)
end

# vQubit
function ITensors.state(::StateName"vId", st::SiteType"vQubit")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
function ITensors.state(::StateName"vX", st::SiteType"vQubit")
    vnames_depwarnmsg()
    return ITensors.state(StateName("X"), st)
end
function ITensors.state(::StateName"vY", st::SiteType"vQubit")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Y"), st)
end
function ITensors.state(::StateName"vZ", st::SiteType"vQubit")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Z"), st)
end
function ITensors.state(::StateName"vH", st::SiteType"vQubit")
    vnames_depwarnmsg()
    return ITensors.state(StateName("H"), st)
end

# vS=1/2
function ITensors.state(::StateName"vSx", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Sx"), st)
end
function ITensors.state(::StateName"vSy", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Sy"), st)
end
function ITensors.state(::StateName"vSz", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Sz"), st)
end
function ITensors.state(::StateName"vX", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("X"), st)
end
function ITensors.state(::StateName"vY", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Y"), st)
end
function ITensors.state(::StateName"vZ", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Z"), st)
end
function ITensors.state(::StateName"vσx", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("σx"), st)
end
function ITensors.state(::StateName"vσy", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("σy"), st)
end
function ITensors.state(::StateName"vσz", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("σz"), st)
end
function ITensors.state(::StateName"vId", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
function ITensors.state(::StateName"vN", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("N"), st)
end
function ITensors.state(::StateName"vecσx", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("σx"), st)
end
function ITensors.state(::StateName"vecσy", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("σy"), st)
end
function ITensors.state(::StateName"vecσz", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("σz"), st)
end
function ITensors.state(::StateName"vecplus", ::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("S+"), st)
end
function ITensors.state(::StateName"vecminus", ::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("S-"), st)
end
function ITensors.state(::StateName"vecN", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("N"), st)
end
function ITensors.state(::StateName"vecId", st::SiteType"vS=1/2")
    vnames_depwarnmsg()
    return ITensors.state(StateName("Id"), st)
end
