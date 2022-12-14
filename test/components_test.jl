using Test
using VirtualInertia
using BlockSystems
import Plots
using OrdinaryDiffEq
import Random
using PlotReferenceTests
set_reference_dir(VirtualInertia)

@testset "ReducedPLL" begin
    rng = Random.MersenneTwister(1)
    pll = Components.ReducedPLL()
    pll = replace_vars(pll; ω_lp=1.32 * 2π*50, Kp=20.0, Ki=2.0)

    p2c = Components.Polar2Cart()
    blk = @connect p2c.(x,y) => pll.(u_r, u_i) outputs=:remaining

    blk = set_input(blk, :mag => 1.0)
    @variables t
    blk = set_input(blk, :arg => 0)

    f = ODEFunction(blk)
    u0 = rand(rng, length(f.syms))
    tspan = (0,10)
    prob = ODEProblem(f,u0,tspan)
    sol = solve(prob, Rodas4())
    Plots.plot(sol)

    # test that angle is close to n*π
    @test abs(rem(sol[end][1], π)) < 0.01
    # test that frequency is near zero
    @test sol[end][2] < 1e-4

    # more complicated test
    pll = Components.ReducedPLL()
    pll = replace_vars(pll; ω_lp=1.32 * 2π*50, Ki=20.0, Kp=3.0)
    p2c = Components.Polar2Cart()
    blk = @connect p2c.(x,y) => pll.(u_r, u_i) outputs=:remaining
    blk = set_input(blk, :mag => 1.0)
    @variables t
    blk = set_input(blk, :arg => 0.1*sin(0.5*t))

    f = ODEFunction(blk)
    u0 = zeros(length(f.syms))
    f.syms .=> u0
    tspan = (0,50)
    prob = ODEProblem(f,u0,tspan)
    sol = solve(prob, Rodas4())

    Plots.plot(sol.t, sol[:ω_pll]; label="tracked freq")
    @reftest "ReducedPLL_tracking_1" Plots.plot!(t->0.1*cos(0.5*t)*0.5; label="input freq")

    Plots.plot(sol.t, sol[:δ_pll]; label="tracked arg")
    @reftest "ReducedPLL_tracking_2" Plots.plot!(t->0.1*sin(0.5*t); label="input arg")
end

@testset "KauraPLL" begin
    rng = Random.MersenneTwister(1)
    pll = Components.KauraPLL()
    pll = replace_vars(pll; ω_lp=500, Kp=20, Ki=3)

    p2c = Components.Polar2Cart()
    blk = @connect p2c.(x,y) => pll.(u_r, u_i) outputs=:remaining

    blk = set_input(blk, :mag => 1.0)
    @variables t
    blk = set_input(blk, :arg => 0)

    f = ODEFunction(blk)
    u0 = rand(rng, length(f.syms))
    u0[5] = 0
    tspan = (0,10)
    prob = ODEProblem(f,u0,tspan)
    sol = solve(prob, Rodas4())
    Plots.plot(sol)

    # test that angle is close to n*π
    @test abs(rem(sol[end][1], π)) < 0.01
    # test that frequency is near zero
    @test sol[end][2] < 1e-4

    # more complicated test
    pll = Components.KauraPLL()
    pll = replace_vars(pll; ω_lp=500, Kp=20., Ki=2)
    p2c = Components.Polar2Cart()
    blk = @connect p2c.(x,y) => pll.(u_r, u_i) outputs=:remaining
    blk = set_input(blk, :mag => 1.0)
    @variables t
    blk = set_input(blk, :arg => 0.1*sin(0.5*t))

    f = ODEFunction(blk)
    u0 = zeros(length(f.syms))
    u0[3] = 1.0
    f.syms .=> u0
    tspan = (0,50)
    prob = ODEProblem(f,u0,tspan)
    sol = solve(prob, Rodas4())

    Plots.plot(sol.t, sol[:ω_pll]; label="tracked freq")
    @reftest "KauraPLL_tracking_1" Plots.plot!(t->0.1*cos(0.5*t)*0.5; label="input freq")

    Plots.plot(sol.t, sol[:δ_pll]; label="tracked arg")
    @reftest "KauraPLL_tracking_2" Plots.plot!(t->0.1*sin(0.5*t); label="input arg")
end

@testset "PT1CurrentSource" begin
    cs = Components.PT1CurrentSource()

    p2c = Components.Polar2Cart()
    blk = @connect p2c.(x,y) => cs.(u_r, u_i) outputs=:remaining

    blk = replace_vars(blk; τ=0.01, P=1, mag=1, Q=0)

    @variables t
    blk = set_input(blk, :arg => 0.5*(t>1))

    f = ODEFunction(blk)
    u0 = zeros(length(f.syms))
    tspan = (0,2)
    prob = ODEProblem(f,u0,tspan)
    sol = solve(prob, Rodas4())
    @reftest "PT1_current_source" Plots.plot(sol)
end

@testset "PerfectCurrentSource" begin
    cs = Components.PerfectCurrentSource()

    gen = generate_io_function(cs; f_inputs=[:u_r, :u_i], f_states=[:i_r, :i_i], f_params=[:P, :Q]);

    rng = Random.MersenneTwister(1)

    for i in 1:100
        u = randn(rng, 2)
        S = randn(rng, 2)
        i = gen.f_oop(u, S, nothing)
        Scomp = Complex(u...) * conj(Complex(i...))

        δP = S[1] .- real(Scomp)
        δQ = S[2] .- imag(Scomp)
        @test δP < 1e-10
        @test δQ < 1e-10
    end
end
