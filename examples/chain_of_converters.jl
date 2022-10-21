using VirtualInertia
using BlockSystems
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
import Plots
using Unitful
using Memoize
using ProgressMeter
#=
# Experiment on virtual inertia
=#

## This parameter set was used in a small EMT experiment with
ω0    = 2π*50u"rad/s"
Sbase = 5u"kW"
Vbase = 230u"V"
Ibase = Sbase/(Vbase) |> u"A"
Cbase = Ibase/Vbase
Lbase = Vbase/Ibase
Rbase = (Vbase^2)/Sbase

## values from 10km
Rline = 0.354u"Ω" / Rbase       |> u"pu"
Lline = 350e-6u"H" / Lbase      |> u"s"
Cline = 2*(12e-6)u"F" / Cbase   |> u"s"

NODES = 5

#=
Topology
```
                           convverters
                  load     3   4   N-1
  second ctrl 1 o---o------o---o ⋯ -o
```
=#

@memoize SerializeDict function solvesystem(;τ_P=0.9e-6, τ_Q=0.9e-6, K_P=1.0, K_Q=1.0, τ_sec=1, tmax=5, τ_load=0.001, n=NODES)
    ## lossless RMS Pi Model Line
    rmsedge = RMSPiLine(R=0, L=ustrip(Lline), C1=ustrip(Cline/2), C2=ustrip(Cline/2))

    ## conv = ODEVertex(VirtualInertia.PT1Source(;τ=0.001),
    conv = ODEVertex(VirtualInertia.PerfectSource(),
                     VirtualInertia.DroopControl(;Q_ref=0,
                                         V_ref=1, ω_ref=0,
                                         τ_P, K_P,
                                         τ_Q, K_Q))

    ## second_ctrl = SecondaryControlCS_PI(; EMT=false, Ki=Ki_sec, Kp=Kp_sec)
    second_ctrl = SecondaryControlCS_PT1(; EMT=false, τ_cs=τ_load, τ=τ_sec)
    second_ctrl.f.params

    load = PT1PLoad(;τ=τ_load)

    ## make sure that we defined all parameters but P_ref for all of the nodes
    load.f.params # P_ref
    conv.f.params # P_ref
    second_ctrl.f.params # P_ref

    g = path_graph(n)

    Nconv = n-2

    nd = network_dynamics([second_ctrl, load, [conv for i in 1:Nconv]...], rmsedge, g)

    uguess = u0guess(nd)
    nd.syms .=> uguess

    pref = [0, -1.0, [1/Nconv for i in 1:Nconv]...]
    p = (pref, nothing)
    ssprob = SteadyStateProblem(nd, uguess, p)
    u0 = solve(ssprob, DynamicSS(Rodas4()))

    nd.syms .=> round.(u0; digits=2)

    function affect(integrator)
        pref_new = copy(pref)
        pref_new[1] += .1 # step up secondary control
        pref_new[2] -= .1 # step down load
        integrator.p = (pref_new, nothing)
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback(0.1, affect)

    tspan = (0.0, tmax)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4(), dtmax=0.01)
end

@nospecialize
plotsym(sol::ODESolution, sym, nodes=1:NODES; mean=false) = plotsym!(plot(), sol, sym, nodes; mean)
plotsym!(sol::ODESolution, sym, nodes=1:NODES; mean=false) = plotsym!(Plots.current(), sol, sym, nodes; mean)
function plotsym!(p::Plots.Plot, sol::ODESolution, sym, nodes=1:NODES; mean=false)
    function names(i)
        if i==1
            " @ secondary control"
        elseif i==2
            " @ load"
        else
            " @ conv$(i-2)"
        end
    end
    if mean
        Plots.plot!(p, meanseries(sol, nodes, sym); label=string(sym)*" mean")
    else
        for i in nodes
            try
                Plots.plot!(p, timeseries(sol, i, sym); label=string(sym)*names(i))
            catch
                @warn "Could not create timeseries for $sym on node $i. Skipped!"
            end
        end
    end
    p
end

function powerloss(sol)
    t, P = timeseries(sol, 1, :Pmeas)
    tidx = findfirst(x->x ≥ 0.1, t)
    t = t[tidx:end]
    P = P[tidx:end]
    Punder = 0.1 .- P
    sum(diff(t) .* Punder[1:end-1])
end
@specialize

## Toggle warnings in block systems
BlockSystems.WARN[] = false

# First we start with with very small values of τ_P and τ_Q (i.e. low invertia)
sol = solvesystem(τ_P=1e-6, τ_Q=1e-6, K_P=1.0, K_Q=1.0, τ_sec=0.8, tmax=5);
@time timeseries(sol,1,:Pmeas);
@time plotsym(sol, :Pmeas)

# lets have a look at different values for the droop term K

VirtualInertia.set_ts_dtmax(0.0005)
plts = Any[]
storage = Vector{Float64}[]
ploss = Float64[]
τ = 0.06
@time for K in [0.5, 1.0, 1.5, 2.0]
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, τ_sec=1, tmax=4);
    p1 = plotsym(sol, :Pmeas, [3:NODES...])
    ## Plots.ylims!(0.95,1.45)
    p2 = plotsym(sol, :imag, 3:NODES)
    p3 = plotsym(sol, :ωmeas, [3:NODES...,1])
    Plots.xlims!(p3, 0.09, 0.2)
    Plots.ylims!(p3, -0.12, 0.01)
    title!(p1, "K = $K")
    push!(plts, Plots.plot(p1, p2, p3, layout=(1,:)))
    push!(storage, needed_storage(sol, 3:NODES))
    push!(ploss, powerloss(sol))
end
plot(plts..., layout=(:,1), size=(1600,1000))

storage
ploss
((sum.(storage) .- ploss)./ploss)*100

#=
For a fixed K, look at different τ
=#
plts = Any[]
storage = Vector{Float64}[]
ploss = Float64[]
K = 1.0
for τ in [1.0, 0.1, 0.01, 0.001]
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, τ_sec=1, tmax=4);
    p1 = plotsym(sol, :Pmeas, [3:NODES...])
    ## Plots.ylims!(0.95,1.45)
    p2 = plotsym(sol, :imag, [3:NODES...])
    ## Plots.ylims!(0.7,0.95)
    p3 = plotsym(sol, :ωmeas, [3:NODES...])
    ## Plots.xlims!(p3, 0.09, 0.2)
    title!(p1, "τ = $τ")
    push!(plts, Plots.plot(p1, p2, p3, layout=(1,:)))
    push!(storage, needed_storage(sol, 3:NODES))
    push!(ploss, powerloss(sol))
end
plot(plts..., layout=(:,1), size=(1600,800))
storage
ploss
((sum.(storage) .- ploss)./ploss)*100

Plots.xlims!(0,1)

#=
Keep rocof constant?
=#
plts = Any[]
@time for K in [0.5, 1.0, 1.5, 2.0]
    innerplts = Any[]
    for M in [0.1, 0.2, 0.3]
        τ = M*K
        sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, τ_sec=1, tmax=2);
        ## p = plotsym(sol, :rocof, [3:NODES...])

        p = plotsym(sol, :ωmeas, [3:NODES...])
        t, ωmean = timeseries(sol, 3, :ωmeas)
        for i in 4:NODES
            _, ω = timeseries(sol, i, :ωmeas)
            ωmean += ω
        end
        ωmean = ωmean / (NODES-2)
        Plots.plot!(t, ωmean, width=5)

        title!(p,"K=$K  M=$M  τ=$τ")
        push!(innerplts, p)
    end
    push!(plts, Plots.plot(innerplts..., layout=(1,:)))
end
plot(plts..., layout=(:,1), size=(1600,1000))

#=
Surface plot
=#
using GLMakie: Makie

Mrange = collect(LinRange(0.1, 0.3, 5))
Krange = collect(LinRange(0.5, 2.0, 5))
rocofs = Array{Float64}(undef, length(Mrange), length(Krange))
nadirs = similar(rocofs)

@time @showprogress for (x, M) in enumerate(Mrange)
    for (y, K) in enumerate(Krange)
        τ = M*K
        sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, τ_sec=1, tmax=2);
        ## p = plotsym(sol, :rocof, [3:NODES...])
        rocofs[x,y] = ROCOF(sol, 3)
        nadirs[x,y] = NADIR(sol, 3)
    end
end

Makie.surface(Mrange, Krange, rocofs, axis=(type=Axis3,ylabel="Droop K", xlabel="virtual mass", zlabel="rocof"))
Makie.surface(Mrange, Krange, nadirs, axis=(type=Axis3,ylabel="Droop K", xlabel="virtual mass", zlabel="nadir"))


#=
surface plots are no good
=#
using DataFrames

results = DataFrame()

Mrange = collect(LinRange(0.1, 0.3, 10))
Krange = collect(LinRange(0.5, 2.0, 5))
for K in Krange
    ## K = 0.5
    ## rocofs = Float64[]
    ## rocoflow = similar(rocofs)
    ## rocofhi = similar(rocofs)
    ## nadirs = similar(rocofs)
    ## nadirlow = similar(rocofs)
    ## nadirhi = similar(rocofs)
    @showprogress for M in Mrange
        τ = M*K
        sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, τ_sec=1, tmax=2);
        ## p = plotsym(sol, :rocof, [3:NODES...])
        roc = ROCOF(sol, 3:NODES)
        nad = NADIR(sol, 3:NODES)
        push!(rocofs, roc[3])
        push!(nadirs, nad[3])
        min, max = extrema(values(roc))
        push!(rocoflow, abs(min - roc[3]))
        push!(rocofhi, abs(max - roc[3]))

        min, max = extrema(values(nad))
        push!(nadirlow, abs(min - nad[3]))
        push!(nadirhi, abs(max - nad[3]))
    end
    Plots.plot!(p1, Mrange, rocofs)#; ribbon=(rocoflow, rocofhi))
    Plots.plot!(p2, Mrange, nadirs; ribbon=(nadirlow, nadirhi))
end


p1 = Plots.plot(; xlabel="M", ylabel="rocof")
p2 = Plots.plot(; xlabel="M", ylabel="nadir")


using DataFrames

df = DataFrame()
Mrange = collect(LinRange(0.1, 0.3, 10))
Krange = collect(LinRange(0.5, 2.0, 5))
@showprogress for K in Krange, M in Mrange
    τ = M*K
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, τ_sec=1, tmax=2);
    ## p = plotsym(sol, :rocof, [3:NODES...])
    for N in 3:NODES
        rocof = ROCOF(sol, N)
        nadir = NADIR(sol, N)
        imax = maximum(timeseries(sol, N, :imag)[2])
        push!(df, (; K, M, τ, N, rocof, nadir, imax))
    end
end

gdf = groupby(df, [:K, :M, :τ])

cmb = combine(gdf,
              AsTable([:rocof, :N]) => (x->DataFrame(x)[x.N .== 3, :rocof]) => :rocof,
              AsTable([:imax, :N]) => (x->DataFrame(x)[x.N .== 3, :imax]) => :imax,
              AsTable([:nadir, :N]) => (x->DataFrame(x)[x.N .== 3, :nadir]) => :nadir,
              :rocof => (x->minimum(abs.(x))) => :rocof_min,
              :rocof => (x->maximum(abs.(x))) => :rocof_max,
              :nadir => (x->minimum(abs.(x))) => :nadir_min,
              :nadir => (x->maximum(abs.(x))) => :nadir_max,
              :imax => (x->minimum(abs.(x))) => :imax_min,
              :imax => (x->maximum(abs.(x))) => :imax_max)

p1 = Plots.plot(title="ROCOF over M")
p2 = Plots.plot(title="nadir over M")
p3 = Plots.plot(title="imax over M")
for i in sort(unique(cmb[!,:K]))
    dfl = cmb[cmb[!,:K] .== i, :]
    Plots.plot!(p1, dfl.M, dfl.rocof; label="K = $i")
    ## ribbon=(dfl.rocof_min-abs.(dfl.rocof), dfl.rocof_max-abs.(dfl.rocof)))
    Plots.plot!(p2, dfl.M, dfl.nadir; label="K = $i")
    Plots.plot!(p3, dfl.M, dfl.imax,;label="K = $i")
end
plot(p1,p2,p3; layout=(:,1), size=(1000,1000))
