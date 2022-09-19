using Test
using EMTSim
using BlockSystems
using Plots
using OrdinaryDiffEq
using SteadyStateDiffEq
using NetworkDynamics
using Graphs
using DiffEqCallbacks


@testset "SecondaryControl_CS" begin
    rmsedge = RMSPiLine(R=0,
                        L=3.30812854442344e-5,
                        C1=0.00012696,
                        C2=0.00012696)
    second_ctrl = SecondaryControlCS(; EMT=false, Ki=5, Kp=0)
    ####
    #### First on slack
    ####
    slack = Slack()

    g = complete_graph(2)
    nd = network_dynamics([second_ctrl, slack], rmsedge, g)
    uguess = u0guess(nd)
    p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 10.0)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4())
    plot(sol)

    ####
    #### Then on droop controller
    ####
    droopPT1 = ODEVertex(EMTSim.PT1Source(τ=0.001),
                         EMTSim.DroopControl(Q_ref=0,
                                             V_ref=1, ω_ref=0,
                                             τ_P=0.01, K_P=0.1,
                                             τ_Q=0.01, K_Q=0.1))

    g = complete_graph(2)
    nd = network_dynamics([second_ctrl, droopPT1], rmsedge, g)

    uguess = u0guess(nd)

    p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 30.0)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4())
    plot(sol)

    plot(timeseries(sol, 1, :Vmag), label="Vmag @ secondary control")
    plot!(timeseries(sol, 2, :Vmag), label="Vmag @ droop")

    plot(timeseries(sol, 1, :Varg), label="Varg @ secondary control")
    plot!(timeseries(sol, 2, :Varg), label="Varg @ droop")

    u0 = sol[end]

    ####
    #### Change injected power and watch pi action
    ####
    function affect(integrator)
        integrator.p = ([-1.0, 1.1], nothing)
        auto_dt_reset!(integrator)
    end
    tspan = (0, 10)
    cb = PresetTimeCallback(1.0, affect)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4())

    plot(timeseries(sol, 1, :Vmag), label="Vmag @ secondary control")
    plot!(timeseries(sol, 2, :Vmag), label="Vmag @ droop")

    plot(timeseries(sol, 1, :Varg), label="Varg @ secondary control")
    plot!(timeseries(sol, 2, :Varg), label="Varg @ droop")

    plot(timeseries(sol, 1, :Pmeas), label="Pmeas @ secondary control")
    plot!(timeseries(sol, 2, :Pmeas), label="Pmeas @ secondary control")

    plot(timeseries(sol, 1, :ω_pll), label="ω_pll @ secondary control")
    plot(twinx(),timeseries(sol, 1, :P_ref_pi), label="Pref @ secondary control")



end
