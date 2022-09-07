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

#=
## RMS Simulation: Droop connected to load
=#

load = ConstPLoad()

# By ommiting the `P_ref=` kw argument the inner Parameter `P_Ref` becomes a
# parameter of the ODEVertex. We may inspect that by calling:
load.f.params


#=
We make P_ref an parameter for the droop too.
=#
droopPT1 = ODEVertex(EMTSim.PT1Source(τ=0.001),
                     EMTSim.DroopControl(Q_ref=0,
                                         V_ref=1, ω_ref=0,
                                         τ_P=0.01, K_P=0.1,
                                         τ_Q=0.01, K_Q=0.1))
droopPT1.f.params

#=
Now we may create the RMS line, ND object, initial guess and find the
fixpoint of the system.
=#
rmsedge = RMSPiLine(R=0, L=ustrip(Lline), C1=ustrip(Cline/2), C2=ustrip(Cline/2))
g = complete_graph(2)
nd = network_dynamics([load, droopPT1], rmsedge, g)
uguess = u0guess(nd)
p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
ssprob = SteadyStateProblem(nd, uguess, p)
u0 = solve(ssprob, DynamicSS(Rodas4()))

#=
To test the inverer response, we need to add a callback to the system which increases
the load at a certain point.
=#
function affect(integrator)
    integrator.p = ([-1.1, 1.0], nothing)
    auto_dt_reset!(integrator)
end
cb = PresetTimeCallback(0.1, affect)

tspan = (0.0, 1.0)
prob = ODEProblem(nd, u0, tspan, p; callback=cb)
sol = solve(prob, Rodas4(), dtmax=0.1)

pmeas1 = plot(timeseries(sol,1,:P_meas); label="P_meas at load")
pmeas2 = plot(timeseries(sol,2,:P_meas); label="P_meas at conv")
plot!(pmeas2, timeseries(sol,2,:P_fil); label="P_fil at conv")

qmeas2 = plot(timeseries(sol,2,:Q_meas); label="Q_meas at conv")
plot!(qmeas2, timeseries(sol,2,:Q_fil); label="Q_fil at conv")

ωplot = plot(timeseries(sol,2,:ω); label="ω at conv")
vabc = plot(timeseries(sol,2,:Va); label="Va at conv")
plot!(timeseries(sol,2,:Vb); label="Vb")
plot!(timeseries(sol,2,:Vc); label="Vc")

Vplot = plot(timeseries(sol,2,:Vmag); label="Vmag at conv")
plot!(Vplot, timeseries(sol,2,:Vmag_ref); label="Vmag setpoint")

allp = [vabc, pmeas1, pmeas2, qmeas2, Vplot, ωplot]
plot(allp..., layout=(6,1), size=(1000,1500))
xlims!(0.099, 0.16)


plot(timeseries(sol,2,:imag))

V = plot(timeseries(sol,1,:Vmag); label="Vmag at load")
plot!(twinx(),timeseries(sol,2,:Vmag); label="Vmag at conv")

Varg= plot(timeseries(sol,1,:Varg); label="Varg at load")
plot!(twinx(), timeseries(sol,2,:Varg); label="Varg at conv")

plot(ωplot, Varg, layout=(2,1))

plot(timeseries(sol,2,:δ))

i = plot(timeseries(sol,1,:imag); label="imag at load")
plot!(i, timeseries(sol,2,:imag); label="imag at conv")
iarg= plot(timeseries(sol,1,:iarg); label="iarg at load")
plot(twinx(), timeseries(sol,2,:imag); label="iarg at conv")


####
#### Now with load
####

load = PT1PLoadEMT(ω0=ustrip(u"rad/s", ω0),
                   C=ustrip(u"s", Cline)/2,
                   τ=1/(2π*50))

droopPT1 = ODEVertex(EMTSim.PT1Source(τ=0.001),
                     EMTSim.DroopControl(Q_ref=0,
                                         V_ref=1, ω_ref=0,
                                         τ_P=0.01, K_P=0.1,
                                         τ_Q=0.01, K_Q=0.1))

edge = EMTRLLine(R=0,#ustrip(u"pu", Rline),
                 L=ustrip(u"s",Lline),
                 ω0=ustrip(u"rad/s", ω0))

g = complete_graph(2)
nd = network_dynamics([load, droopPT1], edge, g)
uguess = u0guess(nd)
p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
tspan = (0.0, 5)
prob = ODEProblem(nd, uguess, tspan, p)
sol = solve(prob, Rosenbrock23())
u0 = sol[end]

function affect(integrator)
    integrator.p = ([-0.9, 1.0], nothing)
    auto_dt_reset!(integrator)
end
cb = PresetTimeCallback(0.1, affect)

tspan = (0.0, 2.0)
prob = ODEProblem(nd, u0, tspan, p; callback=cb)
sol = solve(prob, Rosenbrock23())

vmag = plot(timeseries(sol,1,:Vmag); label="Vmag at load")
plot!(timeseries(sol,2,:Vmag); label="Vmag at Conv")

varg = plot(timeseries(sol,1,:Varg); label="Varg at load")
plot!(timeseries(sol,2,:Varg); label="Varg at Conv")

pmeas1 = plot(timeseries(sol,1,:P_meas); label="P_meas at load")
pmeas2 = plot(timeseries(sol,2,:P_meas); label="P_meas at conv")
plot!(pmeas2, timeseries(sol,2,:P_fil); label="P_fil at conv")

qmeas2 = plot(timeseries(sol,2,:Q_meas); label="Q_meas at conv")
plot!(qmeas2, timeseries(sol,2,:Q_fil); label="Q_fil at conv")

ωplot = plot(timeseries(sol,2,:ω); label="ω at conv")

vabc = plot(timeseries(sol,1,:Va); label="Va at conv")
plot!(timeseries(sol,1,:Vb); label="Vb")
plot!(timeseries(sol,1,:Vc); label="Vc")
xlims!(0.099,0.16)

Vplot = plot(timeseries(sol,2,:Vmag); label="Vmag at conv")
plot!(Vplot, timeseries(sol,2,:Vmag_ref); label="Vmag setpoint")

allp = [vabc, pmeas1, pmeas2, qmeas2, Vplot, ωplot]
plot(allp..., layout=(6,1), size=(1000,1500))
xlims!(0.099, 0.16)

plot!
