# successful demonstration of reconstructing sparse A
# with c_stoner stims!

# big takeaway: Convex.jl is ~200X slower than JuMP.jl!

using Lensman, LinearAlgebra, StatsBase, Gurobi, JuMP
using Convex
import Plots
import Plots: heatmap
m = 1024 # dimension of influence map
n = 32 # number of stimulated neurons
A = randn(m, n)
sparsity = 0.99
base = 2
A[StatsBase.sample(CartesianIndices(A)[:], Int(ceil(n * m * sparsity)), replace=false)] .= 0

function heatmap2(A, B)
    cmin = minimum([A[:]..., B[:]...])
    cmax = maximum([A[:]..., B[:]...])
    p1 = Plots.heatmap(A, clim=(cmin, cmax))
    p2 = Plots.heatmap(B, clim=(cmin, cmax))
    Plots.plot(p1, p2, layout=2)
end

# cstoner approach

stimGroups, groupsPerCell = stonerStimGroups(n, base)
# for sparse matrix A, we try to recover with as few queries as possible
# where A x_1 = y_1 has query x_1 and observation y_1
# X = [ x_1 x_2 ... x_n ]
# Y = [ y_1 y_2 ... y_n ]

nQueries = size(stimGroups, 1)
stimMagnitude = 1

X = zeros(n, nQueries)
for s in 1:nQueries
    neuronsToStim = stimGroups[s, :]
    X[neuronsToStim, s] .= stimMagnitude
end

Y = zeros(m, nQueries)
for s in 1:nQueries
    Y[:, s] = A * X[:, s]
end

## random sampling approach

X_rand = zeros(n, nQueries)
for s in 1:nQueries
    neuronsToStim = StatsBase.sample(1:n, Int(n / 2), replace=false)
    X_rand[neuronsToStim, s] .= stimMagnitude
end

# sum(X_rand, dims=1) # should be n/2 (i.e. 16) for each column

Y_rand = zeros(m, nQueries)
for s in 1:nQueries
    Y_rand[:, s] = A * X_rand[:, s]
end



Â = reconstructA(X, Y, Gurobi.Optimizer);
Â_rand = reconstructA(X_rand, Y_rand, Gurobi.Optimizer);

@show sum(abs.(A .- Â))
@show sum(abs.(A .- Â_rand))

# p = heatmap2(X,Y)
# Plots.savefig(p, "/tmp/XY.png")

