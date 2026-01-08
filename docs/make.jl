using TERRA
using Documenter

DocMeta.setdocmeta!(TERRA, :DocTestSetup, :(using TERRA); recursive=true)

makedocs(;
    modules=[TERRA],
    authors="Amin Taziny",
    sitename="TERRA.jl",
    format=Documenter.HTML(;
        canonical="https://amta3208.github.io/TERRA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amta3208/TERRA.jl",
    devbranch="main",
)
