using VirtualInertia
using Test
using SafeTestsets
using BlockSystems

BlockSystems.WARN[] = false

@safetestset "Components Tests" begin include("components_test.jl") end
@safetestset "Models Tests" begin include("models_test.jl") end
@safetestset "SerializeDict Tests" begin include("SerializeDict_test.jl") end
@safetestset "Solution inspection" begin include("solution_inspection_test.jl") end
@safetestset "Utils" begin include("utils_test.jl") end

BlockSystems.WARN[] = true
