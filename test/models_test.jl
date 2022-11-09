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

using PlotReferenceTests
set_reference_dir(VirtualInertia)

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
    load = ODEVertex(ConstLoad(Q_load=0), [:P_load])
    droopPT1 = ODEVertex(VirtualInertia.PT1Source(τ=0.1),
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

    tspan = (0.0, 0.5)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4(), dtmax=0.0001)

    pmeas1 = Plots.plot(timeseries(sol,1,:_P))
    pmeas2 = Plots.plot(timeseries(sol,2,:_P))
    # power raises to 1.1, filtered slower than actual
    @reftest "droop_P" Plots.plot!(pmeas2, timeseries(sol,2,:P_fil))

    # omega sinks to Δω = K*ΔP = 0.1*0.1 = -0.01
    @reftest "droop_omega" ωplot = Plots.plot(timeseries(sol,2,:ω))

    vabc = Plots.plot(timeseries(sol,2,:_u_a); label="Va at conv")
    Plots.plot!(timeseries(sol,2,:_u_b); label="Vb")
    Plots.plot!(timeseries(sol,2,:_u_c); label="Vc")
    Plots.xlims!(0.098,0.16)

    Vplot = Plots.plot(timeseries(sol,2,:_u_mag))
    Plots.plot!(Vplot, timeseries(sol,2,:PT1Src₊u_ref_mag))
    Plots.plot!(Vplot, timeseries(sol,1,:_u_mag))
    @reftest "droop_voltage_magnitude" Vplot

    argplot = Plots.plot(timeseries(sol,2,:_u_arg); label="uarg at conv")
    Plots.plot!(argplot, timeseries(sol,2,:PT1Src₊u_ref_arg); label="uarg setpoint")
    Plots.plot!(argplot, timeseries(sol,1,:_u_arg); label="uarg @load")
    @reftest "droop_voltage_angle" argplot
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

    Plots.plot(timeseries(sol, 1, :_u_mag))
    Plots.plot!(timeseries(sol, 2, :_u_mag))

    Plots.plot(timeseries(sol, 1, :_u_arg))
    Plots.plot!(timeseries(sol, 2, :_u_arg))

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

    Plots.plot(timeseries(sol, 1, :_u_mag))
    Plots.plot!(timeseries(sol, 2, :_u_mag))

    Plots.plot(timeseries(sol, 1, :_u_arg))
    Plots.plot!(timeseries(sol, 2, :_u_arg))

    Plots.plot(timeseries(sol, 1, :_P))
    @reftest "SecCtrlPI_P" Plots.plot!(timeseries(sol, 2, :_P))

    Plots.plot(timeseries(sol, 1, :ω_pll))
    @reftest "SecCtrlPI_PLL" Plots.plot!(Plots.twinx(),timeseries(sol, 1, :P_ref_pi))
end

@testset "compare RMSPiLine and BSPiLine" begin
    import Random
    Random.seed!(1)
    L = rand()
    R = rand()
    C1 = rand()
    C2 = rand()
    line1 = RMSPiLine(;L, R, C1, C2)
    # Y = G + im B => X = -im/B => B = 2π C
    # Z = R + im X => X = 2π L
    line2 = StaticEdge(BSPiLine(; R, X=2π*50*L, B_src=2π*50*C1, B_dst=2π*50*C2))
    @test line1.sym == line2.sym

    v_src = rand(2)
    v_dst = rand(2)
    e1 = zeros(4)
    e2 = zeros(4)
    line1.f(e1, v_src, v_dst, nothing, 0.0)
    line2.f(e2, v_src, v_dst, nothing, 0.0)
    @test e1 ≈ e2
end

@testset "PT1 and const Load" begin
    slack = ODEVertex(Slack())
    constload = ODEVertex(ConstLoad(), [:P_load, :Q_load])
    pt1load   = ODEVertex(PT1Load(τ=0.001), [:P_load, :Q_load])

    X = ustrip(u"pu", Lline*ω0)
    B = ustrip(u"pu", Cline*ω0)
    line = BSPiLine(;R=0, X, B_src=B, B_dst=B) |> StaticEdge

    g = SimpleGraph(3)
    add_edge!(g,1,2)
    add_edge!(g,1,3)

    nd = network_dynamics([slack,constload,pt1load], line, g)
    u0 = u0guess(nd)
    p = ([(NaN,NaN), (-1, -0.1), (-1,-0.1)] ,nothing)
    prob = ODEProblem(nd, u0, (0,0.1), p)
    sol = solve(prob, Rodas4())

    Plots.plot(timeseries(sol, 1, :_u_mag))
    Plots.plot!(timeseries(sol, 2, :_u_mag))
    @reftest "loads_on_slack_u_mag" Plots.plot!(timeseries(sol, 3, :_u_mag))

    Plots.plot(timeseries(sol, 1, :_u_arg))
    Plots.plot!(timeseries(sol, 2, :_u_arg))
    @reftest "loads_on_slack_u_arg" Plots.plot!(timeseries(sol, 3, :_u_arg))

    Plots.plot(timeseries(sol, 1, :_P))
    Plots.plot!(timeseries(sol, 2, :_P))
    @reftest "loads_on_slack_P" Plots.plot!(timeseries(sol, 3, :_P))

    Plots.plot(timeseries(sol, 1, :_Q))
    Plots.plot!(timeseries(sol, 2, :_Q))
    @reftest "loads_on_slack_Q" Plots.plot!(timeseries(sol, 3, :_Q))
end

@testset "PT1 and const Load (EMT)" begin
    slack = ODEVertex(Slack())
    Cbus = 0.002
    constload = ODEVertex(ConstLoad(EMT=true, C=Cbus, ω0=2π*50), [:P_load, :Q_load])
    pt1load   = ODEVertex(PT1Load(EMT=true, C=Cbus, ω0=2π*50, τ=0.001), [:P_load, :Q_load])

    X = ustrip(u"pu", Lline*ω0)
    B = ustrip(u"pu", Cline*ω0)
    line = BSPiLine(;R=0.0001, X, B_src=B, B_dst=B) |> StaticEdge

    g = SimpleGraph(3)
    add_edge!(g,1,2)
    add_edge!(g,1,3)

    nd = network_dynamics([slack,constload,pt1load], line, g)
    u0 = u0guess(nd)
    p = ([(NaN,NaN), (-1, -0.1), (-1,-0.1)] ,nothing)
    prob = ODEProblem(nd, u0, (0,0.1), p)
    sol = solve(prob, Rodas4(), dtmax=0.0001)

    Plots.plot(timeseries(sol, 1, :_u_mag))
    Plots.plot!(timeseries(sol, 2, :_u_mag))
    @reftest "loads_on_slack_emt_u_mag" Plots.plot!(timeseries(sol, 3, :_u_mag))

    Plots.plot(timeseries(sol, 1, :_u_arg))
    Plots.plot!(timeseries(sol, 2, :_u_arg))
    @reftest "loads_on_slack_emt_u_arg" Plots.plot!(timeseries(sol, 3, :_u_arg))

    Plots.plot(timeseries(sol, 1, :_P))
    Plots.plot!(timeseries(sol, 2, :_P))
    @reftest "loads_on_slack_emt_P" Plots.plot!(timeseries(sol, 3, :_P))

    Plots.plot(timeseries(sol, 1, :_Q))
    Plots.plot!(timeseries(sol, 2, :_Q))
    @reftest "loads_on_slack_emt_Q" Plots.plot!(timeseries(sol, 3, :_Q))

    Plots.plot(timeseries(sol, 2, :_u_a); xlims=(0,0.02))
    Plots.plot!(timeseries(sol, 2, :_u_b))
    @reftest "loads_on_slack_emt_u_abc1" Plots.plot!(timeseries(sol, 2, :_u_c))

    Plots.plot(timeseries(sol, 3, :_u_a); xlims=(0,0.08))
    Plots.plot!(timeseries(sol, 3, :_u_b))
    @reftest "loads_on_slack_emt_u_abc2" Plots.plot!(timeseries(sol, 3, :_u_c))
end

@testset "EMT: Droop on EMT PT1 load" begin
    Cbus = 0.002 # needs to be relative high othervise voltage drop during initialization
    load1 = PT1Load(EMT=true, ω0=ustrip(u"rad/s", ω0),
                   C=Cbus,
                   τ=1/(2π*50)*10e-4,
                   Q_load=0)
    load2 = ConstLoad(EMT=true, ω0=ustrip(u"rad/s", ω0),
                   C=Cbus, Q_load=0)
    load1 = ODEVertex(load1, [:P_load])
    load2 = ODEVertex(load2, [:P_load])

    droopPT1 = ODEVertex(VirtualInertia.PT1Source(;τ=0.001),
                         VirtualInertia.DroopControl(Q_ref=0,
                                             V_ref=1, ω_ref=0,
                                             τ_P=0.01, K_P=0.1,
                                             τ_Q=0.01, K_Q=0.1))

    edge = EMTRLLine(R=ustrip(u"pu", Rline),
                     L=ustrip(u"s",Lline),
                     ω0=ustrip(u"rad/s", ω0)) |> ODEEdge

    g = SimpleGraph(3); add_edge!(g,1,2); add_edge!(g,1,3)
    nd = network_dynamics([droopPT1,load1,load2], edge, g)
    uguess = u0guess(nd)
    p = ([1.0, -.5, -.5], nothing) # node parameters and (empty) line parameters
    tspan = (0.0, 5)
    prob = ODEProblem(nd, uguess, tspan, p)
    sol = solve(prob, Rodas4());
    Plots.plot(sol)
    u0 = sol[end]

    # change the power
    function affect(integrator)
        integrator.p = ([1.0, -.55, -.55], nothing)
        auto_dt_reset!(integrator)
    end
    cb = PresetTimeCallback(0.1, affect)

    tspan = (0.0, 0.5)
    prob = ODEProblem(nd, u0, tspan, p; callback=cb)
    sol = solve(prob, Rodas4(); dtmax=0.0001);

    Plots.plot(timeseries(sol, 1, :_ω))
    Plots.plot!(timeseries(sol, 2, :_ω))
    Plots.plot!(timeseries(sol, 3, :_ω))
    Plots.ylims!(-0.02,0.02)
    @reftest "droop_loads_emt_omega" Plots.xlims!(0.07,0.2)

    Plots.plot(timeseries(sol,1,:_u_mag))
    Plots.plot!(timeseries(sol,2,:_u_mag))
    @reftest "droop_loads_emt_u_mag" Plots.plot!(timeseries(sol,3,:_u_mag))

    Plots.plot(timeseries(sol,1,:_u_arg))
    Plots.plot!(timeseries(sol,2,:_u_arg))
    @reftest "droop_loads_emt_u_arg" Plots.plot!(timeseries(sol,3,:_u_arg))
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

    Plots.plot(timeseries(sol, 1, :_P))
    @reftest "SecCtrlPT1_P" Plots.plot!(timeseries(sol, 2, :_P))
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

    @reftest "Syncrhonverter_omega" Plots.plot(timeseries(sol, 2, :ω))
    @reftest "Syncrhonverter_mfif" Plots.plot(timeseries(sol, 2, :MfIf))

    @reftest "Syncrhonverter_P" Plots.plot(timeseries(sol, 2, :_P))
    @reftest "Syncrhonverter_Q" Plots.plot!(timeseries(sol, 2, :_Q))

    @reftest "Syncrhonverter_u_mag" Plots.plot(timeseries(sol, 2, :_u_mag))

    @reftest "Syncrhonverter_u_arg" Plots.plot(timeseries(sol, 2, :_u_arg))
end

@testset "Synchronverter on load with change" begin
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
    prob = ODEProblem(nd, u0, tspan, p, callback=cb)
    sol = solve(prob, Rodas4(), dtmax=0.1)

    Plots.plot(timeseries(sol,1,:_P))
    @reftest "Synchronverter_loadchange_P" Plots.plot!(timeseries(sol,2,:_P))

    # Plots.plot(timeseries(sol,2,:MfIf))
    # Plots.plot(timeseries(sol,2,:Te))
    # Plots.plot(timeseries(sol,2,:Q))

    @reftest "Synchronverter_loadchange_omega" ωplot = Plots.plot(timeseries(sol,2,:ω))
end
