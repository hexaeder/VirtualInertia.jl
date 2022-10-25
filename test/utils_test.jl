using Test
using OrdinaryDiffEq
using DiffEqCallbacks
using VirtualInertia

@testset "PRecord test" begin
    tspan = (0, 1)
    f(du,u,p,t) = du = p[1][1] * u
    params = (rand(4), ["foo", "bar"])

    prob = ODEProblem(f, [1], tspan, params)
    precord = VirtualInertia.PRecord(prob)

    function affect(integrator)
        newp = deepcopy(integrator.p)
        newp[1][1] = integrator.t
        newp[2][1] = "bas at $(integrator.t)"
        integrator.p = (newp[1], newp[2])
        record!(precord, integrator)
    end
    cb1 = PresetTimeCallback(0.1, affect)
    cb2 = PresetTimeCallback(0.5, affect)

    sol = solve(prob, Tsit5(); callback=CallbackSet(cb1, cb2));

    @test precord(0) == precord(0.05)
    @test precord(0) != precord(0.1)
    @test precord(0.1) == precord(0.4)
    @test precord(0.1) != precord(1)

    @test precord(0.1)[1][1] == precord(0.1; direction=:right)[1][1] != precord(0.1; direction=:left)[1][1]

    precord([0,0.1,0.1,0.5,0.5,1])

end
