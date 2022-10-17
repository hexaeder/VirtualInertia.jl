using CSV
using DataFrames
using Graphs
using MetaGraphs
using Missings
using Unitful
using GeometryBasics
using NetworkLayout

export import_system, describe_nodes, describe_edges, bustype, is_static_state

const DATA_DIR = abspath(@__DIR__, "..", "data")

"""
    bustype(i::Int)
    bustype(s::String)
    bustype(s::Symbol)

Helper to convert the different representations of bus types into eachother.
Base representation is Symbol:

    :load   == bustype(1) == bustype("PQ") == bustype(:PQ)
    :gen    == bustype(2) == bustype("PV") == bustype(:PV)
    :syncon == bustype(3) == bustype("Ref") == bustype(:Ref)
"""
bustype(i::Int) = bustype(Val(i))
bustype(s) = bustype(Val(Symbol(s)))
bustype(::Union{Val{:PQ},  Val{1}, Val{:load}})   = :load
bustype(::Union{Val{:PV},  Val{2}, Val{:gen}})    = :gen
bustype(::Union{Val{:Ref}, Val{3}, Val{:syncon}}) = :syncon

"""
    import_system(sym::Symbol; kwargs...)::MetaGraph

Main entry point to load the systems. New systems should overload this function. Known implementations
- `import_system(:rtsgmlc)`: loads the GMLC update for the rts96

Those functions return a `MetaGraph` which has properties attached to the Nodes/Edges.

Graph properties:
- `Pbase` : Base power for PU
optional:
- `NodeProps` : tuple of node property names which should appear first in describe functions
- `EdgeProps` : tuple of edge property names which should appear first in describe functions

Node properties:
- `P` : active power in PU
- `Q` : reactive power in PU
- `Vbase` : base voltage in kV
- `Vm` : voltage magnitude in PU
optional:
- `inertia` : inertia in MJ/MW (see H in Wikipedia page)
- `damping` : damping factor γ for swing equations
- `timeconstant` : time cosntant `τ` for dynamic loads
- `x`, `y` : position of Bus for plotting purposes
- `Va` : voltage angle in rad

Edge properties:
- `rating` : short term emergency rating in PU
- `R` : resistance in PU
- `X` : reactance in PU
"""
import_system(sym::Symbol; kwargs...) = import_system(Val(sym); kwargs...)

"""
    describe_nodes(g::MetaGraph; firstcols=Vector{String}())

Returns DataFrame with all the node meta data.
"""
function describe_nodes(g::MetaGraph; firstcols=Vector{String}())
    df = DataFrame(; n=Int[])
    for n in 1:nv(g)
        row = push!(props(g, n), :n=>n)
        push!(df, row, cols=:union)
    end
    firstcols = String.(firstcols) # convert to string if given as Symbol
    has_prop(g, :NodeProps) &&
        append!(firstcols, String.(get_prop(g, :NodeProps)))
    append!(firstcols, names(df)) |> unique!
    select!(df, firstcols)
end

"""
    describe_nodes(g::MetaGraph; firstcols=Vector{String}())

Returns DataFrame with all the edge meta data.
"""
function describe_edges(g::MetaGraph; firstcols=Vector{String}())
    df = DataFrame(; src=Int[], dst=Int[])
    for e in edges(g)
        row = push!(props(g, e), :src=>e.src, :dst=>e.dst)
        push!(df, row, cols=:union)
    end
    firstcols = String.(firstcols) # convert to string if given as Symbol
    has_prop(g, :EdgeProps) &&
        append!(firstcols, String.(get_prop(g, :EdgeProps)))
    append!(firstcols, names(df)) |> unique!
    select!(df, firstcols)
end

import MetaGraphs: set_prop!, get_prop, has_prop

KEY_ITER = Union{AbstractUnitRange,Vector,AbstractEdgeIter}
"""
    set_prop!(g, keys::Iterable, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for differet identifiers `keys`.
"""
function set_prop!(g, keys::KEY_ITER, prop::Symbol, vals::Vector)
    length(keys) == length(vals) || throw(ArgumentError("keys and vals needs to be of same length!"))
    for (k, val) in zip(keys, vals)
        if !ismissing(val)
            set_prop!(g, k, prop, val)
        end
    end
end

set_prop!(g, keys::KEY_ITER, p::Symbol, val) = set_prop!(g, keys, p, [val for v in keys])

"""
    get_prop(g, keys::Iterable, prop::Symbol)

Get same property `prop` with different values `vals` for differet keys.
"""
function get_prop(g, keys::KEY_ITER, prop::Symbol)
    [has_prop(g, k, prop) ? get_prop(g, k, prop) : missing for k in keys]
end

has_prop(g, keys::KEY_ITER, prop::Symbol) = all(k -> has_prop(g, k, prop), keys)

"""
    import_system(:rtsgmlc)

Import the RTS-GMLC system as a MetaGraph.
"""
function import_system(::Val{:rtsgmlc}; losses=false, kwargs...)
    @info "Import system RTS-GMLC"
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data    = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    gen_data    = CSV.read(joinpath(data,"gen.csv"), DataFrame)
    branch_data = CSV.read(joinpath(data,"branch.csv"), DataFrame)

    sort!(bus_data, "Bus ID")
    N = nrow(bus_data)
    g = MetaGraph(N)

    baseP = 100u"MW"
    set_prop!(g, :Pbase, baseP)
    set_prop!(g, :NodeProps, [:n, :type, :P_load, :Q_load, :P_inj, :Q_inj,:H, :id,  :Vm, :Vbase,
                              :Va ])
    set_prop!(g, :EdgeProps, [:src, :dst, :R, :X, :rating, :id])

    set_prop!(g, 1:N, :id, bus_data."Bus ID")
    set_prop!(g, 1:N, :Vm, bus_data."V Mag"u"pu")
    set_prop!(g, 1:N, :Va, bus_data."V Angle"u"°")
    set_prop!(g, 1:N, :Vbase, bus_data."BaseKV"u"kV")
    # set_prop!(g, 1:N, :x, bus_data."lng")
    # set_prop!(g, 1:N, :y, bus_data."lat")
    # bustypes = bustype.(bus_data."Bus Type")
    # set_prop!(g, 1:N, :type, bustypes)

    pload = - bus_data."MW Load"u"MW" ./ baseP * u"pu"
    qload = - bus_data."MVAR Load"u"MW" ./ baseP * u"pu"
    set_prop!(g, 1:N, :P_load, pload)
    set_prop!(g, 1:N, :Q_load, qload)

    for (n, busid) in enumerate(bus_data."Bus ID")
        P_load = get_prop(g, n, :P_load)
        Q_load = get_prop(g, n, :Q_load)
        # select attached generators
        controltype = bus_data."Bus Type"[n]
        generators = gen_data[gen_data."Bus ID".==busid, ["Category", "MW Inj", "MVAR Inj", "Inertia MJ/MW"]]
        if isempty(generators)
            P_inj = 0.0u"pu"
            Q_inj = 0.0u"pu"
        else # has generators
            P_inj = sum(generators."MW Inj")u"MW" / baseP * u"pu"
            Q_inj = sum(generators."MVAR Inj")u"MW" / baseP * u"pu"
            set_prop!(g, n, :P_inj, P_inj)
            set_prop!(g, n, :Q_inj, Q_inj)
        end

        if controltype ∈ ("PV", "Ref")
            if bus_data."Bus ID"[n] ∈ [114, 214, 314]
                type = :syncon
                set_prop!(g, n, :H, 5u"MJ/MW")
            else
                type = :gen
                inertia = sum(generators."Inertia MJ/MW")u"MJ/MW"
                @assert !iszero(P_inj) "Generator $n doese not inject power?"
                @assert !iszero(inertia) "Generator $n doese not have inrtia?"
                set_prop!(g, n, :H, inertia)
            end
        elseif controltype == "PQ"
            type = :load
        else
            error("Found unkonown control type $controltype")
        end

        set_prop!(g, n, :type, type)
        # set_prop!(g, n, :P, P_inj - P_load)
        # set_prop!(g, n, :Q, Q_inj - Q_load)
    end

    # set thepe intertia for sync condenser
    # scidxs = findall(s -> !ismissing(s) && iszero(s), describe_nodes(g).H)
    # @assert bus_data."Bus ID"[scidxs] == [114, 214, 314]
    # set_prop!(g, scidxs, :H, 5u"MJ/MW")

    for row in eachrow(branch_data)
        src = findfirst(x->x==row."From Bus", bus_data."Bus ID")
        dst = findfirst(x->x==row."To Bus", bus_data."Bus ID")
        propertys = Dict(:R => row."R"u"pu",
                         :X => row."X"u"pu",
                         :B => row."B"u"pu",
                         :rating => row."STE Rating"u"MW"/baseP * u"pu",
                         :id => row."UID")
        add_edge!(g, src, dst, propertys)
    end

    set_gmlc_pos_relaxed!(g)

    losses || lossless!(g)
    set_missing_props!(g; kwargs...)

    return g
end

export balance_power!, lossless!
function balance_power!(network)
    nodes = describe_nodes(network)
    imbalance = sum(skipmissing(nodes.P_inj)) + sum(nodes.P_load)
    if imbalance ≈ 0
        println("already balanced!")
        return network
    end
    genidx = findall(!ismissing, nodes.P_inj)

    relative_inj = nodes.P_inj[genidx] ./ sum(nodes.P_inj[genidx])

    newp = copy(nodes.P_inj)
    newp[genidx] .-= relative_inj .* imbalance

    pdiff = sum(skipmissing(newp)) + sum(nodes.P_load)
    @assert isapprox(pdiff, 0, atol=1e-8) "Could not balance power! Sum is $pdiff)"

    set_prop!(network, 1:nv(network), :P_inj, newp)
end

function lossless!(network)
    balance_power!(network)
    set_prop!(network, edges(network), :R, zeros(ne(network)) * u"pu")
end

function set_missing_props!(network; damping = nothing, tconst = nothing)
    nodes = describe_nodes(network)
    if !isnothing(damping)
        idxs = findall(x -> x === :gen || x === :syncon, bustype.(nodes.type))
        set_prop!(network, idxs, :damping, damping)
    end
    if !isnothing(tconst)
        idxs = findall(x -> x === :load, bustype.(nodes.type))
        set_prop!(network, idxs, :timeconst, tconst)
    end
end

function set_gmlc_pos_relaxed!(g)
    # set location property based on rts gmlc data
    data = joinpath(DATA_DIR, "RTS-GMLC")
    bus_data    = CSV.read(joinpath(data,"bus.csv"), DataFrame)
    x = bus_data."lng"
    y = bus_data."lat"
    xn = 10*(x .- minimum(x))./(maximum(x)-minimum(x)).-5
    yn = 10*(y .- minimum(y))./(maximum(y)-minimum(y)).-5
    pos2 = spring(g, initialpos=Point2f0.(zip(xn, yn)))

    pos2[21] += Point2(-.3,.25)
    pos2[15] += Point2(-.5,.0)
    pos2[17] += Point2(-.3,-.3)
    pos2[24] += Point2(-.3,-.3)
    pos2[18] += Point2(-.2,-.0)

    set_prop!(g, 1:nv(g), :pos, pos2)
end
