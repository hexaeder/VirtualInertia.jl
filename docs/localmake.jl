#! julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
Pkg.instantiate()
Pkg.update()

using LiveServer
@async serve(dir=joinpath(@__DIR__, "build"))

run = true
while run
    try
        include("make.jl")
    catch e
        @warn "There was an error during `make.jl`" e
    end

    println("Run again? Enter! Exit with 'q'.")
    if readline() == "q"
        global run = false
    end
end
