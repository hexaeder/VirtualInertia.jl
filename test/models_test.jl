using Test
using VirtualInertia
using BlockSystems
import Plots
using OrdinaryDiffEq
using SteadyStateDiffEq
using NetworkDynamics
using Graphs
using DiffEqCallbacks
using Unitful

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

@testset "Doop Control test on step load" begin
    load = ODEVertex(ConstPLoad(), [:P_ref])
    droopPT1 = ODEVertex(VirtualInertia.PT1Source(τ=0.001),
                         VirtualInertia.DroopControl(Q_ref=0,
                                             V_ref=1, ω_ref=0,
                                             τ_P=0.01, K_P=0.1,
                                             τ_Q=0.01, K_Q=0.1))
    @test droopPT1.f.params == [:P_ref]

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


    VirtualInertia.set_ts_dtmax(0.0005)
    pmeas1 = Plots.plot(timeseries(sol,1,:P_meas); label="P_meas at load")
    pmeas2 = Plots.plot(timeseries(sol,2,:P_meas); label="P_meas at conv")
    Plots.plot!(pmeas2, timeseries(sol,2,:P_fil); label="P_fil at conv")

    qmeas2 = Plots.plot(timeseries(sol,2,:Q_meas); label="Q_meas at conv")
    Plots.plot!(qmeas2, timeseries(sol,2,:Q_fil); label="Q_fil at conv")

    ωplot = Plots.plot(timeseries(sol,2,:ω); label="ω at conv")
    vabc = Plots.plot(timeseries(sol,2,:Va); label="Va at conv")
    Plots.plot!(timeseries(sol,2,:Vb); label="Vb")
    Plots.plot!(timeseries(sol,2,:Vc); label="Vc")

    Vplot = Plots.plot(timeseries(sol,2,:Vmag); label="Vmag at conv")
    Plots.plot!(Vplot, timeseries(sol,2,:Vmag_ref); label="Vmag setpoint")

    allp = [vabc, pmeas1, pmeas2, qmeas2, Vplot, ωplot]
    Plots.plot(allp..., layout=(6,1), size=(1000,1500))
    Plots.xlims!(0.099, 0.16)
end

@testset "SecondaryControlCS_PI" begin
    rmsedge = RMSPiLine(R=0,
                        L=3.30812854442344e-5,
                        C1=0.00012696,
                        C2=0.00012696)
    second_ctrl = SecondaryControlCS_PI(; EMT=false, Ki=5, Kp=0)
    second_ctrl = ODEVertex(second_ctrl, [:P_ref])

    ####
    #### First on slack
    ####
    slack = ODEVertex(Slack())

    g = complete_graph(2)
    nd = network_dynamics([second_ctrl, slack], rmsedge, g)
    uguess = u0guess(nd)
    p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 10.0)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4())
    Plots.plot(sol)

    ####
    #### Then on droop controller
    ####
    droopPT1 = ODEVertex(VirtualInertia.PT1Source(τ=0.001),
                         VirtualInertia.DroopControl(Q_ref=0,
                                             V_ref=1, ω_ref=0,
                                             τ_P=0.01, K_P=0.1,
                                             τ_Q=0.01, K_Q=0.1));

    g = complete_graph(2)
    nd = network_dynamics([second_ctrl, droopPT1], rmsedge, g)

    uguess = u0guess(nd)

    p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 30.0)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4())
    Plots.plot(sol)

    Plots.plot(timeseries(sol, 1, :Vmag), label="Vmag @ secondary control")
    Plots.plot!(timeseries(sol, 2, :Vmag), label="Vmag @ droop")

    Plots.plot(timeseries(sol, 1, :Varg), label="Varg @ secondary control")
    Plots.plot!(timeseries(sol, 2, :Varg), label="Varg @ droop")

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

    Plots.plot(timeseries(sol, 1, :Vmag), label="Vmag @ secondary control")
    Plots.plot!(timeseries(sol, 2, :Vmag), label="Vmag @ droop")

    Plots.plot(timeseries(sol, 1, :Varg), label="Varg @ secondary control")
    Plots.plot!(timeseries(sol, 2, :Varg), label="Varg @ droop")

    Plots.plot(timeseries(sol, 1, :Pmeas), label="Pmeas @ secondary control")
    Plots.plot!(timeseries(sol, 2, :Pmeas), label="Pmeas @ secondary control")

    Plots.plot(timeseries(sol, 1, :ω_pll), label="ω_pll @ secondary control")
    Plots.plot!(Plots.twinx(),timeseries(sol, 1, :P_ref_pi), label="Pref @ secondary control")
end

@testset "PT1 and const Load" begin
    slack = ODEVertex(Slack())
    # load = PT1Load(EMT=true, τ=0.01, C=1e-5, ω0=2π*50)
    # load = ConstLoad(EMT=true, C=1e-5, ω0=2π*50)
    load = ConstLoad()
    load = PT1Load(τ=0.001)
    loadv = ODEVertex(load,[:P_load, :Q_load])

    # line = RMSPiLine(R=0, L=3.3e-5, C1=1e-4, C2=1e-4)
    # line = RMSPiLine(R=0, L=0.8/(2π*50), C1=1e-4/(2π*50), C2=1e-4/(2π*50))
    line = BSPiLine(;R=0, X=0.8, B_src=0.05, B_dst=0.05) |> StaticEdge
    # line = BSPiLine(;R=0, X=3.3e-5*(2π*50), B_src=1e-4*(2π*50), B_dst=1e-4*(2π*50)) |> StaticEdge

    g = SimpleGraph(2)
    add_edge!(g,1,2)
    nd = network_dynamics([slack, loadv], line, g)
    u0 = u0guess(nd)
    p = ([(NaN,NaN), (-1.0, 0)] ,nothing)
    prob = ODEProblem(nd, u0, (0,0.1), p)
    @time sol = solve(prob, Rodas4());

    Plots.plot(sol)
end

@testset "EMT: Droop on EMT PT1 load" begin
    load = PT1Load(EMT=true, ω0=ustrip(u"rad/s", ω0),
                   C=ustrip(u"s", Cline)/2,
                   τ=1/(2π*50),
                   Q_load=0)
    load = ODEVertex(load, [:P_load])

    droopPT1 = ODEVertex(VirtualInertia.PT1Source(;τ=0.001),
                         VirtualInertia.DroopControl(Q_ref=0,
                                             V_ref=1, ω_ref=0,
                                             τ_P=0.01, K_P=0.1,
                                             τ_Q=0.01, K_Q=0.1))

    edge = EMTRLLine(R=0,#ustrip(u"pu", Rline),
                     L=ustrip(u"s",Lline),
                     ω0=ustrip(u"rad/s", ω0)) |> ODEEdge

    g = complete_graph(2)
    nd = network_dynamics([load, droopPT1], edge, g)
    uguess = u0guess(nd)
    p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 5)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4())
    Plots.plot(sol)
    u0 = sol[end]

    function affect(integrator)
        integrator.p = ([-0.9, 1.0], nothing)
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback(0.1, affect)

    tspan = (0.0, 2.0)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4())

    vmag = Plots.plot(timeseries(sol,1,:Vmag); label="Vmag at load")
    Plots.plot!(timeseries(sol,2,:Vmag); label="Vmag at Conv")

    varg = Plots.plot(timeseries(sol,1,:Varg); label="Varg at load")
    Plots.plot!(timeseries(sol,2,:Varg); label="Varg at Conv")

    pmeas1 = Plots.plot(timeseries(sol,1,:Pmeas); label="P_meas at load")
    pmeas2 = Plots.plot(timeseries(sol,2,:Pmeas); label="P_meas at conv")
    Plots.plot!(pmeas2, timeseries(sol,2,:P_fil); label="P_fil at conv")

    qmeas2 = Plots.plot(timeseries(sol,2,:Q_meas); label="Q_meas at conv")
    Plots.plot!(qmeas2, timeseries(sol,2,:Q_fil); label="Q_fil at conv")

    ωplot = Plots.plot(timeseries(sol,2,:ω); label="ω at conv")

    VirtualInertia.set_ts_dtmax(0.00001)
    # vabc = Plots.plot(timeseries(sol,1,:Va); label="Va at conv")
    # Plots.plot!(timeseries(sol,1,:Vb); label="Vb")
    # Plots.plot!(timeseries(sol,1,:Vc); label="Vc")
    vabc = Plots.plot(timeseries(sol,2,:ia); label="ia at conv")
    Plots.plot!(timeseries(sol,2,:ib); label="ib")
    Plots.plot!(timeseries(sol,2,:ic); label="ic")
    Plots.xlims!(0.09,0.14)
    # Plots.ylims!(0.7, 0.9)


    Vplot = Plots.plot(timeseries(sol,2,:Vmag); label="Vmag at conv")
    Plots.plot!(Vplot, timeseries(sol,2,:Vmag_ref); label="Vmag setpoint")

    allp = [vabc, pmeas1, pmeas2, qmeas2, Vplot, ωplot]
    Plots.plot(allp..., layout=(6,1), size=(1000,1500))
    Plots.xlims!(0.09, 0.13)
end

@testset "SecondaryControlCS_PT1" begin
    rmsedge = RMSPiLine(R=0,
                        L=3.30812854442344e-5,
                        C1=0.00012696,
                        C2=0.00012696)
    second_ctrl = SecondaryControlCS_PT1(; EMT=false, τ_cs=0.001, τ=1.0)
    second_ctrl = ODEVertex(second_ctrl, [:P_ref])
    second_ctrl.f.params
    ####
    #### First on slack
    ####
    slack = ODEVertex(Slack())

    g = complete_graph(2)
    nd = network_dynamics([second_ctrl, slack], rmsedge, g)
    uguess = u0guess(nd)
    p = ([-1.0, 1.0], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 10.0)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4())
    Plots.plot(sol)

    Plots.plot(timeseries(sol, 2, :Pmeas))
end

@testset "Synchronverter on slack" begin
    slack = ODEVertex(Slack());
    syncv = ODEVertex(VirtualInertia.PT1Source(τ=0.001),
                      VirtualInertia.Synchronverter(P_ref=1, Q_ref=0.1,
                                            V_ref=1,
                                            J=0.02, Dp=0.1, Dq=10.1, Kv=1,
                                            ω0=ustrip(ω0)));

    rmsedge = RMSPiLine(R=0, L=ustrip(Lline), C1=ustrip(Cline/2), C2=ustrip(Cline/2))
    g = complete_graph(2)
    nd = network_dynamics([slack, syncv], rmsedge, g)
    uguess = u0guess(nd)
    prob = ODEProblem(nd, uguess, (0, 5))
    sol = solve(prob, Rodas4());

    Plots.plot(timeseries(sol, 2, :ω)...)
    Plots.plot(timeseries(sol, 2, :MfIf)...)

    Plots.plot(timeseries(sol, 2, :Pmeas)...)
    Plots.plot(timeseries(sol, 2, :Qmeas)...)

    Plots.plot(timeseries(sol, 2, :Vmag)...)
    Plots.plot(timeseries(sol, 2, :Varg)...)

    Plots.plot(timeseries(sol, 2, :Vmag_ref)...)
    Plots.plot(timeseries(sol, 2, :Varg_ref)...)
end

@testset "Synchronverter on load with chagne" begin
    # load = ConstPLoad()
    load = ODEVertex(PT1Load(τ=0.001, Q_load=0), [:P_load])
    @test load.f.params == [:P_load]

    syncv = ODEVertex(VirtualInertia.PT1Source(τ=0.001),
                      VirtualInertia.Synchronverter(Q_ref=0,
                                            V_ref=1,
                                            J=0.02, Dp=0.1, Dq=1.0, Kv=1,
                                            ω0=ustrip(ω0)));
    @test syncv.f.params == [:P_ref]

    rmsedge = RMSPiLine(R=0, L=ustrip(Lline), C1=ustrip(Cline/2), C2=ustrip(Cline/2))

    g = complete_graph(2)
    nd = network_dynamics([load, syncv], rmsedge, g)
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

    tspan = (0.0, 2.)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4(), dtmax=0.1)

    Plots.plot(timeseries(sol,1,:Pmeas); label="P_meas at load")
    Plots.plot!(timeseries(sol,2,:Pmeas); label="P_meas at conv")

    Plots.plot(timeseries(sol,2,:MfIf))
    Plots.plot(timeseries(sol,2,:Te))
    Plots.plot(timeseries(sol,2,:Q))

    ωplot = Plots.plot(timeseries(sol,2,:ω); label="ω at conv")
end
