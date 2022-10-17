using VirtualInertia
using Test
using SafeTestsets
using BlockSystems

BlockSystems.WARN[] = false

@safetestset "Components Tests" begin include("components_test.jl") end
@safetestset "Models Tests" begin include("models_test.jl") end
@safetestset "SerializeDict Tests" begin include("SerializeDict_test.jl") end

BlockSystems.WARN[] = true
