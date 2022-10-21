using VirtualInertia
using MetaGraphs
using Graphs
using BlockSystems
using Unitful
using NetworkDynamics
using OrdinaryDiffEq
using SteadyStateDiffEq
using DiffEqCallbacks
# import Plots
using GLMakie: Makie
using DataFrames
import Base.Iterators

#####
##### Import Network
#####
network = import_system(:rtsgmlc)
nodes_df = describe_nodes(network)
edges_df = describe_edges(network)

#####
##### Line Model
#####
lineblk = BSPiLine(; B_src=:B, B_dst=:B)
staticpiline = StaticEdge(lineblk, [:R, :X, :B])

#####
##### Inverter Model
#####
inner = PerfectSource()
outer = DroopControl(;V_ref=1, ω_ref=0, τ_P=:τ, τ_Q=:τ, K_P=:K, K_Q=:K)
# here i've renamed K_Q and K_P in a way that they are the same

invblk = Inverter(inner, outer; name=:inv)
# we can add a load to the voltage source
invloadblk = VSwithLoad(invblk)

inverter = ODEVertex(invloadblk, [:P_load, :Q_load, :P_ref, :Q_ref, :τ, :K])

#####
##### Load Model
#####
loadblk = PT1Load(τ=0.001)
# loadblk = ConstLoad()
# loadblk = PT1Load(EMT=true, τ=0.001, C=1e-4, ω0=2π*50)
# loadblk = ConstLoad(EMT=true, C=1e-4, ω0=2π*50)
load = ODEVertex(loadblk, [:P_load, :Q_load])

#####
##### Create ND.jl parameter construct
#####
# line parameters are (R, L, C)
R = ustrip.(u"pu", edges_df.R)
X = ustrip.(u"pu", edges_df.X)
B = ustrip.(u"pu", edges_df.B) ./ 2
edgeP = zip(R, X, B) |> collect

# node parameters are
# Inverters: Pload, Qload, Pref, Qref, τ, K
gens = filter(n -> n.type ∈ (:gen, :syncon), nodes_df)

τinv = 0.01
Kinv = 0.1

gen_para = collect(zip(ustrip.(u"pu", gens.P_load),
                       ustrip.(u"pu", gens.Q_load),
                       ustrip.(u"pu", gens.P_inj),
                       ustrip.(u"pu", gens.Q_inj),
                       Iterators.cycle(τinv),
                       Iterators.cycle(Kinv)))

# loadsPload, Qload for loads
loads = filter(n -> n.type == :load, nodes_df)
load_para = collect(zip(ustrip.(u"pu", loads.P_load),
                        ustrip.(u"pu", loads.Q_load),
                        Iterators.cycle(NaN),
                        Iterators.cycle(NaN),
                        Iterators.cycle(NaN),
                        Iterators.cycle(NaN)))

nodeP = [Tuple(NaN for i in 1:6) for i in 1:length(nodes_df.n)]
nodeP[gens.n]  = gen_para
nodeP[loads.n] = load_para

# FINALLY create the parameter struct
params = (nodeP, edgeP)

#####
##### Crate NetworkDynamics model
#####
vertices = Vector{ODEVertex}(undef, length(nodes_df.n))
vertices[gens.n] .= Ref(inverter);
vertices[loads.n] .= Ref(load);

# solve the steady state problem
nd = network_dynamics(vertices, staticpiline, SimpleGraph(network))
uguess = u0guess(nd)
ssprob = SteadyStateProblem(nd, uguess, params)
u0 = solve(ssprob, DynamicSS(Rodas4()))

#####
##### Disturb system
#####
largest = findmax(x -> ismissing(x) ? 0 : x, nodes_df.P_inj)[2]

disturbed = 10
disturbance = -0.1

function affect(integrator)
    newp = deepcopy(integrator.p)
    distp    = newp[1][disturbed] |> collect
    distp[1] += disturbance   # increase load
    largestp = newp[1][largest] |> collect
    largestp[3] -= disturbance # add this to p inj
    newp[1][disturbed] = Tuple(distp)
    newp[1][largest] = Tuple(largestp)

    integrator.p = (newp[1], newp[2])
    auto_dt_reset!(integrator)
end
cb = PresetTimeCallback(0.1, affect)

tspan = (0, 100)
prob = ODEProblem(nd, u0, tspan, params; callback=cb)
sol = solve(prob, Rodas4());

plotsym(sol, :ω, largest)
plotsym(sol, :P_fil, largest)
plotsym(sol, :Varg, largest)
plotsym(sol, :δ, largest)

Makie.plot(timeseries(sol, 1, :ω))
#TODO why is the frequenccy not zero but the angle stabilizes?!

Makie.plot(1:10000,y->sin(y))
