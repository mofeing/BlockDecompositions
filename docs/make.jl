using BlockDecompositions
using Documenter

DocMeta.setdocmeta!(BlockDecompositions, :DocTestSetup, :(using BlockDecompositions); recursive=true)

makedocs(;
    modules=[BlockDecompositions],
    authors="Sergio Sánchez Ramírez <sergio.sanchez.ramirez@bsc.es> and contributors",
    repo="https://github.com/Sergio Sánchez Ramírez/BlockDecompositions.jl/blob/{commit}{path}#{line}",
    sitename="BlockDecompositions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
