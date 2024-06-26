using LindbladVectorizedTensors
using ITensors
using ITensorMPS
using Documenter

DocMeta.setdocmeta!(
    LindbladVectorizedTensors,
    :DocTestSetup,
    :(using LindbladVectorizedTensors);
    recursive=true,
)

makedocs(;
    modules=[LindbladVectorizedTensors],
    checkdocs=:exported,
    authors="Davide Ferracin <davide.ferracin@protonmail.com> and contributors",
    sitename="LindbladVectorizedTensors.jl",
    format=Documenter.HTML(;
        canonical="https://phaerrax.github.io/LindbladVectorizedTensors.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "New site types" => "included_site_types.md",
        "Examples" => "gksl_example.md",
        "Reference" => "physical_operators.md",
    ],
)

deploydocs(;
    branch="gh-pages",
    repo="github.com/phaerrax/LindbladVectorizedTensors.jl",
    devbranch="main",
)
