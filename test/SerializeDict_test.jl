using VirtualInertia
using Test

@testset "test SerializeDict" begin
    dir = tempname(pwd())
    sd = SerializeDict(dir)
    empty!(sd)

    sd[1] = "foo"
    sd[2] = "bar"

    @test sd[1] == "foo"
    @test sd[2] == "bar"
    @test_throws KeyError sd[3]

    sd[4] = "bax"

    @test length(sd) == 3
    delete!(sd, 4)
    @test length(sd) == 2
    delete!(sd, 41)
    @test length(sd) == 2

    @test get(sd, 1, "baz") == "foo"
    @test get(sd, 3, "baz") == "baz"
    @test_throws KeyError sd[3]
    @test get!(sd, 3, "baz") == "baz"
    @test sd[3] == "baz"

    @test get!(sd, 7) do
        "foo"*"bar"
    end == "foobar"

    @test sd[7] == "foobar"
    rm(dir; recursive=true)
end
