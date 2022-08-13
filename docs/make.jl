using CompressedSparseArrays
using Documenter

DocMeta.setdocmeta!(CompressedSparseArrays, :DocTestSetup, :(using CompressedSparseArrays); recursive=true)

makedocs(;
    modules=[CompressedSparseArrays],
    authors="Will Kimmerer <kimmerer@mit.edu> and contributors",
    repo="https://github.com/Wimmerer/CompressedSparseArrays.jl/blob/{commit}{path}#{line}",
    sitename="CompressedSparseArrays.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Wimmerer.github.io/CompressedSparseArrays.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Wimmerer/CompressedSparseArrays.jl",
    devbranch="main",
)
