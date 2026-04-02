function vecHv_depwarnmsg(name)
    warnmsg = "The \"vec$name\" and \"Hv$name\" site types are deprecated and will be \
    removed in a future version of this package. Please use the \"v$name\" site type \
    instead."
    Base.depwarn(warnmsg, :space; force=true)
end

ITensors.alias(::SiteType"HvS=1/2") = SiteType"vS=1/2"()
ITensors.alias(::SiteType"vecS=1/2") = SiteType"vS=1/2"()

function ITensors.space(st::SiteType"HvS=1/2")
    vecHv_depwarnmsg("S=1/2")
    return ITensors.space(ITensors.alias(st))
end
function ITensors.space(st::SiteType"vecS=1/2")
    vecHv_depwarnmsg("S=1/2")
    return ITensors.space(ITensors.alias(st))
end

function ITensors.val(vn::ValName, st::SiteType"HvS=1/2")
    vecHv_depwarnmsg("S=1/2")
    return ITensors.val(vn, ITensors.alias(st))
end
function ITensors.val(vn::ValName, st::SiteType"vecS=1/2")
    vecHv_depwarnmsg("S=1/2")
    return ITensors.val(vn, ITensors.alias(st))
end

function ITensors.state(sn::StateName, st::SiteType"vecS=1/2"; kwargs...)
    vecHv_depwarnmsg("S=1/2")
    return ITensors.state(sn, ITensors.alias(st); kwargs...)
end
function ITensors.state(sn::StateName, st::SiteType"HvS=1/2"; kwargs...)
    vecHv_depwarnmsg("S=1/2")
    return ITensors.state(sn, ITensors.alias(st); kwargs...)
end

function ITensors.op(on::OpName, st::SiteType"vecS=1/2"; kwargs...)
    vecHv_depwarnmsg("S=1/2")
    return ITensors.op(on, ITensors.alias(st); kwargs...)
end
function ITensors.op(on::OpName, st::SiteType"HvS=1/2"; kwargs...)
    vecHv_depwarnmsg("S=1/2")
    return ITensors.op(on, ITensors.alias(st); kwargs...)
end
