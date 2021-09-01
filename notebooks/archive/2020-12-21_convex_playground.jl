# successful demonstration of reconstructing sparse A
# with c_stoner stims!

# big takeaway: Convex.jl is ~200X slower than JuMP.jl!

using Convex, Lensman, COSMO, SCS, LinearAlgebra, StatsBase, Gurobi,
    JuMP
import Plots
import Plots: heatmap
m = 1024 # dimension of influence map
n = 32 # number of stimulated neurons
A = randn(m,n)
sparsity = 0.99
base = 2
A[StatsBase.sample(CartesianIndices(A)[:], Int(ceil(n*m*sparsity)), replace=false)] .= 0

function heatmap2(A,B)
    cmin = minimum([A[:]..., B[:]...])
    cmax = maximum([A[:]..., B[:]...])
    p1 = Plots.heatmap(A,clim=(cmin,cmax))
    p2 = Plots.heatmap(B,clim=(cmin,cmax))
    Plots.plot(p1,p2,layout=2)
end

stimGroups, groupsPerCell = stonerStimGroups(n, base)
# for sparse matrix A, we try to recover with as few queries as possible
# where A x_1 = y_1 has query x_1 and observation y_1
# X = [ x_1 x_2 ... x_n ]
# Y = [ y_1 y_2 ... y_n ]

nQueries = size(stimGroups,1)
stimMagnitude = 1

X = zeros(n,nQueries)
for s in 1:nQueries
    neuronsToStim = stimGroups[s,:]
    X[neuronsToStim,s] .= stimMagnitude
end

Y = zeros(m,nQueries)
for s in 1:nQueries
    Y[:,s] = A * X[:,s]
end

# p = heatmap2(X,Y)
# Plots.savefig(p, "/tmp/XY.png")


function benchConvex(X,Y; optimizer=() -> COSMO.Optimizer(verbose = false))
    m,n = size(Y,1), size(X,1)
    M = Convex.Variable(m,n)
    problem = minimize(sum(abs.(M)), M*X == Y)
    Convex.solve!(problem, optimizer)
    evaluate(M)
end

function benchJump(X,Y)
    m,n = size(Y,1), size(X,1)
    model = Model(Gurobi.Optimizer)
    @variable(model, Â[1:m,1:n])
    @variable(model, l1)
    @constraint(model, [l1; vec(Â)] in MOI.NormOneCone(m*n))
    @constraint(model, con, Â*X .== Y)
    @objective(model, Min, l1)
    optimize!(model)
    value.(Â)
end
# all single threaded..?!? (slow to due convex.jl constraints...)
# println("bench COSMO") # 532 seconds
# @time benchConvex(X,Y)
# println("bench SCS") # 264 seconds
# @time benchConvex(X,Y, optimizer=() -> SCS.Optimizer(verbose = false))
println("bench Gurobi") # 144 seconds
@time Ahat = benchConvex(X,Y, optimizer=() -> Gurobi.Optimizer())
@time Â = benchJump(X,Y); # 1.03 seconds!!!
# JuMP.jl is >100x faster than Convex.jl

@time Â = reconstructA(X,Y, optimizer=SCS.Optimizer); # 85s
@time Â = reconstructA(X,Y, optimizer=COSMO.Optimizer); # 15.3s
@time Â = reconstructA(X,Y, optimizer=Gurobi.Optimizer); # 0.88s

# variable performance depending on A--
# sometimes a different solver gets closer..?
# (seems to be random)
@show sum(abs.(A .- Â)) # (1st) 43.418 (2nd) 50.99
@show sum(abs.(A))
@show sum(abs.(A*X .- Â*X))
@show sum(abs.(Ahat))
@show sum(abs.(Â))
@show sum(abs.(A .- Ahat)) # (1st) 39.433 (2nd) 55.60
@show sum(abs.(A*X .- Ahat*X))
heatmap2(A,value.(Â))
##

# Convex.solve!(problem, () -> SCS.Optimizer(max_iters = 10000, verbose = false))

M = Convex.Variable(n,n)
problem = minimize(sum(abs.(M)), M*X == Y)
# Convex.solve!(problem, () -> SCS.Optimizer(max_iters = 10000, verbose = false))
Convex.solve!(problem, () -> COSMO.Optimizer(verbose = false))
@show norm(evaluate(M),1)
@show norm(A,1)
# A*X ≈ evaluate(M)*X
# p = heatmap2(A,evaluate(M))
# Plots.savefig(p, "/tmp/XY.png")
##
pX = Plots.heatmap(X, title="X")
pY = Plots.heatmap(Y, title="Y")
cmin = minimum([A[:]..., Â[:]...])
cmax = maximum([A[:]..., Â[:]...])
pA = Plots.heatmap(A,clim=(cmin,cmax), title="A")
pAhat = Plots.heatmap(Â,clim=(cmin,cmax), title="Ahat")
cmin = minimum([(A*X)[:]..., (Â*X)[:]...])
cmax = maximum([(A*X)[:]..., (Â*X)[:]...])
pAX = Plots.heatmap(A*X,clim=(cmin,cmax), title="A*X")
pAhatX = Plots.heatmap(Â*X,clim=(cmin,cmax), title="Ahat*X")
p = Plots.plot(pX, pA, pAX, pY,  pAhat, pAhatX, layout=(2,3))
##
# Plots.savefig(p, "/tmp/not_optimal.png")
