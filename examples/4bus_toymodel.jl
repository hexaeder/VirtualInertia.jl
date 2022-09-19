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

using EMTSim
using BlockSystems
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using Plots
import CairoMakie
using GraphMakie
using Unitful

#=
Topology
```
  second ctrl 1 o-------o 2 load -2 -> -2.5
                |\      |
                | \     |
                |  \    |
                |   \   |
                |    \  |
                |     \ |
                |      \|
        conv2 4 o-------o 3 conv1 +1
```
=#

function solvesystem(;τ_P=0.9e-6, τ_Q=0.9e-6, K_P=1.0, K_Q=1.0, Ki_sec=1, Kp_sec=0, tmax=5, τ_load=0.001)
    # lossless RMS Pi Model Line
    rmsedge = RMSPiLine(R=0, L=ustrip(Lline), C1=ustrip(Cline/2), C2=ustrip(Cline/2))

    # conv1 = ODEVertex(EMTSim.PT1Source(;τ=0.001),
    conv1 = ODEVertex(EMTSim.PerfectSource(),
                      EMTSim.DroopControl(;Q_ref=0,
                                          V_ref=1, ω_ref=0,
                                          τ_P, K_P,
                                          τ_Q, K_Q))

    # conv2 = ODEVertex(EMTSim.PT1Source(;τ=0.001),
    conv2 = ODEVertex(EMTSim.PerfectSource(),
                      EMTSim.DroopControl(;Q_ref=0,
                                          V_ref=1, ω_ref=0,
                                          τ_P, K_P,
                                          τ_Q, K_Q))

    second_ctrl = SecondaryControlCS(; EMT=false, Ki=Ki_sec, Kp=Kp_sec)

    load = PT1PLoad(;τ=τ_load)

    # make sure that we defined all parameters but P_ref for all of the nodes
    load.f.params # P_ref
    conv1.f.params # P_ref
    conv2.f.params # P_ref
    second_ctrl.f.params # P_ref

    # g = PathGraph(4)
    # add_edge!(g, 4, 1)
    # add_edge!(g, 1, 3)
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 2)

    nd = network_dynamics([second_ctrl, load, conv1, conv2], rmsedge, g)
    uguess = u0guess(nd)
    nd.syms .=> uguess

    pref = [0, -2.0, 1.0, 1.0]
    p = (pref, nothing)
    ssprob = SteadyStateProblem(nd, uguess, p)
    u0 = solve(ssprob, DynamicSS(Rodas4()))

    nd.syms .=> round.(u0; digits=2)

    function affect(integrator)
        pref_new = [0, -2.5, 1.0, 1.0]
        integrator.p = (pref_new, nothing)
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback(0.1, affect)

    tspan = (0.0, tmax)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4())
end

plotsym(sol::ODESolution, sym, nodes=[1,2,3,4]) = plotsym!(plot(), sol, sym, nodes)
plotsym!(sol::ODESolution, sym, nodes=[1,2,3,4]) = plotsym!(Plots.current(), sol, sym, nodes)
function plotsym!(p::Plots.Plot, sol::ODESolution, sym, nodes=[1,2,3,4])
    names = Dict(1 => " @ secondary control",
                 2 => " @ load",
                 3 => " @ conv1",
                 4 => " @ conv2")
    for i in nodes
        try
            plot!(p, timeseries(sol, i, sym); label=string(sym)*names[i])
        catch
            @warn "Could not create timeseries for $sym on node $i. Skipped!"
        end
    end
    p
end

function needed_storage(sol, idx)
    node_p = sol.prob.p[1][idx]
    psyms = EMTSim._getwrapper(sol.prob.f, idx).params
    pidx = findfirst(isequal(:P_ref), psyms)

    t, P = timeseries(sol, idx, :Pmeas)
    Pover = P .- node_p[pidx]
    sum(diff(t) .* Pover[1:end-1])
end

# First we start with with very small values of τ_P and τ_Q (i.e. low invertia)
sol = solvesystem(τ_P=1e-6, τ_Q=1e-6, K_P=1.0, K_Q=1.0, Ki_sec=1, Kp_sec=0, tmax=5);
p = plotsym(sol, :Pmeas, [1,3,4])

# lets have a look at different values for the droop term K
EMTSim.TS_DTMAX 0.0005
EMTSim.set_ts_dtmax(0.0005)
plts = Any[]
for K in [0, 0.1, 0.8, 1.2]
    sol = solvesystem(τ_P=1e-6, τ_Q=1e-6, K_P=K, K_Q=K, Ki_sec=1, Kp_sec=0, tmax=3);
    p1 = plotsym(sol, :Pmeas, [3,4,1])
    # ylims!(0.95,1.45)
    p2 = plotsym(sol, :imag, [3,4])
    p3 = plotsym(sol, :ωmeas, [1,3,4])
    xlims!(p3, 0.09, 0.2)
    title!(p1, "K = $K")
    storage3 = needed_storage(sol, 3)
    storage4 = needed_storage(sol, 3)
    title!(p2, "Conv1 bat=$storage3")
    title!(p3, "Conv2 bat=$storage4")
    push!(plts, plot(p1, p2, p3, layout=(1,:)))
end
plot(plts..., layout=(:,1), size=(1600,1000))

#=
For a fixed K, look at different τ
=#
plts = Any[]
K = 0.8
for τ in [1.3, 0.1, 0.01, 0.001]
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, Ki_sec=1, Kp_sec=0, tmax=12.5);
    p1 = plotsym(sol, :Pmeas, [3,4])
    # ylims!(0.95,1.45)
    p2 = plotsym(sol, :imag, [3,4])
    # ylims!(0.7,0.95)
    p3 = plotsym(sol, :ωmeas, [3,4])
    xlims!(p3, 0.09, 0.2)
    title!(p1, "τ = $τ")
    storage3 = needed_storage(sol, 3)
    storage4 = needed_storage(sol, 4)
    title!(p2, "Conv1 bat=$storage3")
    title!(p3, "Conv2 bat=$storage4")
    push!(plts, plot(p1, p2, p3, layout=(1,:)))
end
plot(plts..., layout=(:,1), size=(1600,800))


#=
Theoretical value for storage needs (what is lost until secondary control kicks in)
=#
int = Any[]
for τ in [1.3, 0.1, 0.01, 0.001]
    K=0.8
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, Ki_sec=1, Kp_sec=0, tmax=10.0);
    t, P = timeseries(sol, 1, :Pmeas)
    Pover = 0.5 .- P
    fidx = findfirst(x->x>0.1, t)
    t = t[fidx:end]
    Pover = Pover[fidx:end]
    push!(int,  sum(diff(t) .* Pover[1:end-1]))
end

p = plot(legend=:bottomright)
for τ in [1.3, 0.1, 0.01, 0.001]
    K=0.8
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, Ki_sec=1, Kp_sec=0, tmax=20.0);
    plotsym!(p, sol, :Pmeas, 1)
end
current()

p = plot(legend=:bottomright)
for τ in [1.3, 0.1, 0.01, 0.001]
    K=0.8
    sol = solvesystem(τ_P=τ, τ_Q=τ, K_P=K, K_Q=K, Ki_sec=1, Kp_sec=0, tmax=20.0);
    plotsym!(p, sol, :ωmeas, 1)
end
current()
