using Test
using VirtualInertia
using NetworkDynamics
using OrdinaryDiffEq
using Graphs

@testset "test _group_states" begin
    using VirtualInertia: _group_states
    states = [:x, :y, :u_r, :u_i, :K]
    _group_states(states)
    using LinearAlgebra
    norm((1,2))
    norm([1,2])
end

@testset "getstate" begin
    load = ODEVertex(VirtualInertia.ConstLoad(Q_load=0,P_load = -1))
    droopPT1 = ODEVertex(VirtualInertia.PT1Source(τ=0.001),
                         VirtualInertia.DroopControl(Q_ref=0, P_ref=1,
                                             V_ref=1, ω_ref=0,
                                             τ_P=0.01, K_P=0.1,
                                             τ_Q=0.01, K_Q=0.1))

    rmsedge = RMSPiLine(R=0, L=0.04, C1=1, C2=1)

    g = complete_graph(2)
    nd = network_dynamics([load, droopPT1], rmsedge, g)
    uguess = u0guess(nd)
    prob = ODEProblem(nd, uguess, (0,1))
    sol = solve(prob, Rodas4())

    blockstates(sol, 1)
    getstate(sol, 0.0,nothing, 1, :u_r)
    getstate(sol, 0.0,nothing, 1, :_u_r)
    getstate(sol, 0.0,nothing, 1, :u_i)
    getstate(sol, 0.0,nothing, 1, :_u_i)
    getstate(sol, 0.0,nothing, 1, :u_mag)
    getstate(sol, 0.0,nothing, 1, :_u_mag)
    getstate(sol, 0.0,nothing, 1, :u_arg)
    getstate(sol, 0.0,nothing, 1, :_u_arg)

    getstate(sol, 0.0,nothing, 1, :_P)
    getstate(sol, 0.0,nothing, 1, :_Q)
    getstate(sol, 0.0,nothing, 1, :_ω)
    getstate(sol, 0.0,nothing, 1, :_rocof)
    getstate(sol, 0.0,nothing, 1, :_u_a)
    getstate(sol, 0.0,nothing, 1, :_i_mag)
end
