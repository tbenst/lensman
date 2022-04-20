using SafeTestsets, Test
println("Time to complete `using Lensman`")
@time using Lensman

@time begin
# @time @safetestset "Algorithm Tests" begin minclude("test_algorithms.jl") end
@time @safetestset "main" begin include("main.jl") end
println("Total test time:")
end