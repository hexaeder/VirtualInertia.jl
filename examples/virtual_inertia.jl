#=
# Experiments on virtual inertia
In this notebook I'd like to explore some first steps regarding the comparison of virtual
inertia in different grid forming inverter strategies.

## Experimental Setup
For now, the setup is quite simple. I connect a single load to the grid.
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


#=
Super sophisticated initial guess.
=#
u0guess(nd::ODEFunction) = u0guess.(nd.syms)
function u0guess(s::Symbol)
    s = string(s)
    if occursin(r"^u_r", s)
        1.0
    elseif occursin(r"^A", s)
        1.0
    else
        0.0
    end
end

using EMTSim
using BlockSystems
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using Plots
using Unitful

slack = Slack()

droopPT1 = ODEVertex(EMTSim.PT1Source(τ=0.01),
                     EMTSim.DroopControl(P_ref=1, Q_ref=0,
                                         V_ref=1, ω_ref=0,
                                         τ_P=0.01, K_P=0.01,
                                         τ_Q=0.01, K_Q=0.01))

rmsedge = RMSPiLine(R=0, L=ustrip(Lline), C1=ustrip(Cline/2), C2=ustrip(Cline/2))

g = complete_graph(2)
nd = network_dynamics([slack, droopPT1], rmsedge, g)
u0 = u0guess(nd)

# @btime nd($(similar(u0)), u0, nothing, 0.0)

tspan = (0, 10)
prob = ODEProblem(nd, u0, tspan)
sol = solve(prob, Rodas4())

plot(timeseries(sol, 2, :Vmag); label="Vmag")
plot!(timeseries(sol, 2, :Vmag_ref); label="Vmag_ref")

plot(timeseries(sol, 2, :Va); label="Va")
plot!(timeseries(sol, 2, :Vb); label="Vb")
plot!(timeseries(sol, 2, :Vc); label="Vc")

plot(timeseries(sol, 2, :Varg); label="Varg")
plot!(timeseries(sol, 2, :Varg_ref); label="Varg_ref")

plot(timeseries(sol, 2, :P); label="P")

plot!(timeseries(sol, 2, :A))
