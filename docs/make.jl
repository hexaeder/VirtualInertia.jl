using EMTSim
using Documenter

DocMeta.setdocmeta!(EMTSim, :DocTestSetup, :(using EMTSim); recursive=true)

makedocs(;
    modules=[EMTSim],
    authors="Hans Würfel <git@wuerfel.io> and contributors",
    repo="https://github.com/hexaeder/EMTSim.jl/blob/{commit}{path}#{line}",
    sitename="EMTSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hexaeder.github.io/EMTSim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hexaeder/EMTSim.jl",
    devbranch="main",
)
