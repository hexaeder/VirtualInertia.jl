#=
# Disconnect load from slack
In this example, we examine the voltage transient of the system
after a single load of constant P is disconnected from a slack.
=#
using EMTSim
using BlockSystems
using ModelingToolkit
using NetworkDynamics
using Graphs
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using Plots
using Unitful
using CSV
using DataFrames

#=
## Constants and unit stuff.

!!! Note
    The grid is modelled in a gloabal dq-Frame without zero component. Since we
    might need to introduce *local* dq frames for the inverters too, the global
    frame is denoted by `_r` and `_i` (for real and imag part) annotations while
    the local frames use `_d` and `_q`.
=#
ω0    = 2π*50u"rad/s"
Sbase = 300u"MW"
Vbase = 110u"kV" #* sqrt(2/3)
Ibase = Sbase/(Vbase)#* √3) # why the √3 ?
Cbase = Ibase/Vbase
Lbase = Vbase/Ibase
Rbase = (Vbase^2)/Sbase

Rline = 1u"Ω" / Rbase           |> u"pu"
Pload = -300u"MW" / Sbase       |> u"pu"
Rload = (1u"pu")^2 / Pload      |> u"pu" # R=U^2/P
Cline = (2e-6)u"F" / Cbase      |> u"s"
Lline = (1/100π)u"H" / Lbase    |> u"s"

nothing#hide
#=
## Slack Bus
The slack bus is modelled as a node with ``\dot(u) = 0``, i'll keep the initial voltage forever.
=#
@variables t u_r(t) u_i(t)
dt = Differential(t)

slackblock = IOBlock([dt(u_r) ~ 0, dt(u_i) ~ 0], [], [u_r, u_i]; name=:slack)

# create ODE Vertex from this block
slack = ODEVertex(slackblock)
nothing#hide

#=
## Constant R Load
I tried to model an constan R load first using the algebraic current equation
```math
i = -\frac 1 R \, u
```
This however, did not match the powerfactory results. The power factory load is
purely ohmic (no phase shift) but constant in P. Since the voltage on bus 2 is slighly below
1 pu, we can not trivially calculate the corresponding resistance.

```
@variables t i_r(t) i_i(t)
@parameters u_r(t) u_i(t) R
loadblock = IOBlock([i_r ~ -1/R * u_r,
                     i_i ~ -1/R * u_i],
                    [u_r, u_i], [i_r, i_i],
                    name=:load)

# The load is used as a current source in a `BusBar`
busblock = BusBar(loadblock; name=:loadbus)
busblock = set_p(busblock, Dict(:C=>ustrip(u"s", Cline), :ω0=>ustrip(u"rad/s", ω0)))

# The `BusBar` block can be used to create an ODE Vertex
load = ODEVertex(busblock, [:load₊R])
```
=#


#=
## Constant P Load
For a constant load P we may use the algebric current equation
```math
i_r = P \frac{u_r}{u_r^2 + u_i^2}\\
i_i = P \frac{u_i}{u_r^2 + u_i^2}
```
However, this leads to instability in the model and is hard to initialize. I guess it is not a
good idea to follow each oszillation in node/condensator voltage!

To circumvent the problem, we model the P load using a lowpass filter

```math
\dot{i_r} = \frac{1}{\tau}\left(P \frac{u_r}{u_r^2 + u_i^2} - i_r\right)\\
\dot{i_i} = \frac{1}{\tau}\left(P \frac{u_i}{u_r^2 + u_i^2} - i_i\right)
```
whose time constant ``\tau = \frac{1}{\omega_0}`` is chosen in a way, that
disturbances higher than the nominal grid frequency are damped.
=#
@variables t i_r(t) i_i(t)
@parameters u_r(t) u_i(t) P τ
loadblock = IOBlock([dt(i_r) ~ ustrip(u"rad/s", ω0)*(P * u_r/(u_r^2 + u_i^2) - i_r),
                     dt(i_i) ~ ustrip(u"rad/s", ω0)*(P * u_i/(u_r^2 + u_i^2) - i_i)],
                    [u_r, u_i], [i_r, i_i],
                    name=:load)

# The load is used as a current source in a `BusBar`
busblock = BusBar(loadblock; name=:loadbus)
busblock = set_p(busblock, Dict(:C=>ustrip(u"s", Cline), :ω0=>ustrip(u"rad/s", ω0)))

# The `BusBar` block can be used to create an ODE Vertex
load = ODEVertex(busblock, [:load₊P])

nothing#hide

#=
## ODE edge for the RL Line
=#
@variables t i_r(t) i_i(t)
@parameters u_r_src(t) u_i_src(t) u_r_dst(t) u_i_dst(t) R L ω

lineblock = IOBlock([dt(i_r) ~  ω * i_i  - R/L * i_r + 1/L*(u_r_src - u_r_dst),
                     dt(i_i) ~ -ω * i_r  - R/L * i_i + 1/L*(u_i_src - u_i_dst)],
                    [u_r_src, u_i_src, u_r_dst, u_i_dst],
                    [i_r, i_i],
                    name=:RLLine)
lineblock = set_p(lineblock, Dict(:R=>NoUnits(Rline), :L=>ustrip(u"s", Lline), :ω=>ustrip(u"rad/s", ω0)))

# This block can be used to create an ODEEdge
edge = ODEEdge(lineblock)

nothing#hide
#=
## Set up the network
=#
g = SimpleGraph(2)
add_edge!(g, 1, 2)
nd = network_dynamics([slack, load], edge, g)

nothing#hide
#=
We use SteadyStateDiffeq to find the steady state based on an initial guess. The initial guess just contains
voltage at both nodes. The load value is given as a parameter to node 2.
=#
uguess = zeros(length(nd.syms))
uguess[[1,3]] .= 1.0 # set the d component of the slack to 1 from the beginning
p = ([0, NoUnits(Pload)], nothing)
ssprob = SteadyStateProblem(nd, uguess, p)
u0 = solve(ssprob, DynamicSS(AutoTsit5(Rosenbrock23())))

nothing#hide
#=
Let's ignore this warning for now. I am not sure why the steady state is not completly steady. it nearly is.

## Simulation
We simulate a disconnection of the load at `t=0.1 s`. This is done using a callback. The callback does 2 things:
it sets the reference point P to zero, it sets the currents flowing through the load to zero.
=#
tspan = (0.0, 0.124)
affect = function (integrator)
    integrator.p = ([0, 0.0], nothing)
    integrator.u[[5,6]] .= 0.0
end
cb = PresetTimeCallback(0.1, affect)
prob = ODEProblem(nd, u0, tspan, p; callback=cb)
sol = solve(prob, Tsit5(),dtmax=0.00001)
nothing #hide

# We may now transformation the results back to a,b,c frame
a,b,c = Tdqinv(sol.t, sol[3,:], sol[4,:])
nothing#hide
#=
## Comparison with power factory results
Read the Power Facory data and plot for reference

``V_{base}`` is RMS phase-phase voltage.


```math
V_{star} = \frac{1}{\sqrt{3}} *  V_{base}\\
\hat{V} = \sqrt{\frac{2}{3}} * V_{base}
```

Therefore, i have to multiply the `a`, `b` and `c` results by `sqrt(3/2)`
=#
df = CSV.read(joinpath(dirname(pathof(EMTSim)), "..", "data","PowerFactory", "Test_EMT.csv"), DataFrame, skipto=3,
              header=[:t, :u_1_a, :u_1_b, :u_1_c, :u_2_a, :u_2_b, :u_2_c])

xlims = (0.099,0.124)
plot(df.t, df.u_2_a; label="PowerFactory A", xlims, color=:lightgray,  linewidth=5)
plot!(sol.t, a.*sqrt(3/2); label="u_2_a", xlims)

plot!(df.t, df.u_2_b; label="PowerFactory B", color=:lightgray, linewidth=5)
plot!(sol.t, b.*sqrt(3/2); label="u_2_b")

plot!(df.t, df.u_2_c; label="PowerFactory C", color=:lightgray, linewidth=5)
plot!(sol.t, c.*sqrt(3/2); label="u_2_c")

#=
In this graph, the PowerFactory solution is light gray and thick behind our solution.

... lets have a closer look
=#

xlims!(0.0995,0.105)
