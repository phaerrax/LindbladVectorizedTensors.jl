# New operators
# -------------
function ITensors.op(on::OpName"X", st::SiteType"Boson", d::Int)
    return 1 / sqrt(2) *
           (ITensors.op(OpName("Adag"), st, d) + ITensors.op(OpName("A"), st, d))
end
function ITensors.op(on::OpName"Y", st::SiteType"Boson", d::Int)
    return im / sqrt(2) *
           (ITensors.op(OpName("Adag"), st, d) - ITensors.op(OpName("A"), st, d))
end
