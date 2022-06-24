using EMTSim
using Documenter
using Literate

DocMeta.setdocmeta!(EMTSim, :DocTestSetup, :(using EMTSim); recursive=true)

# generate examples
example_dir = abspath(joinpath(@__DIR__, "..", "examples"))
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

@info "Create Markdown files from examples"

examples = filter(contains(r".jl$"), readdir(example_dir))
exlinks = Pair{String,String}[]
for example in examples
    path = joinpath(example_dir, example)
    Literate.markdown(path, outdir)

    barename = example[begin:end-3]
    push!(exlinks, barename => joinpath("generated", barename*".md"))
end

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
        "Examples" => exlinks
    ],
)

deploydocs(;
    repo="github.com/hexaeder/EMTSim.jl",
    devbranch="main",
)
