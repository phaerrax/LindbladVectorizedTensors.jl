function trace(x::MPS)
    vid = MPS(siteinds(x), "vId")
    return dot(vid, x)
end

function vectorization_expvalues(; atol=1e-12)
    # Test that the vectorization is done correctly by measuring the expectation values
    # of some observables on a random state.
    sites = siteinds("Qubit", 4)
    x = random_mps(sites; linkdims=4)

    x_vec = vec_purestate_densitymatrix(SiteType("Qubit"), x)
    sites_vec = siteinds(x_vec)

    y_exp = [dot(x, apply(op("Y", sites, i), x)) for i in 1:length(x)]
    y_exp_vec = [trace(apply(op("Y⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]

    h_exp = [dot(x, apply(op("H", sites, i), x)) for i in 1:length(x)]
    h_exp_vec = [trace(apply(op("H⋅", sites_vec, i), x_vec)) for i in 1:length(x_vec)]

    angle = pi * rand()
    cp_exp = [
        dot(x, apply(op("CPhase", sites, i, i + 2; ϕ=angle), x)) for
        i in 1:length(x) if i + 2 <= length(x)
    ]
    cp_exp_vec = [
        trace(apply(op("CPhase⋅", sites_vec, i, i + 2; ϕ=angle), x_vec)) for
        i in 1:length(x_vec) if i + 2 <= length(x_vec)
    ]

    ccx_exp = [
        dot(x, apply(op("CCNOT", sites, i, i + 1, i + 2), x)) for
        i in 1:length(x) if i + 2 <= length(x)
    ]
    ccx_exp_vec = [
        trace(apply(op("CCNOT⋅", sites_vec, i, i + 1, i + 2), x_vec)) for
        i in 1:length(x_vec) if i + 2 <= length(x_vec)
    ]

    return (
        norm(h_exp - h_exp_vec) < atol,
        norm(h_exp - h_exp_vec) < atol,
        norm(cp_exp - cp_exp_vec) < atol,
        norm(ccx_exp - ccx_exp_vec) < atol,
    )
end
