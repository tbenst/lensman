# successful demonstration of reconstructing sparse A
# with c_stoner stims!

using Convex, Lensman, COSMO, SCS, LinearAlgebra
import Plots
n = 32
A = zeros(n,n)
base=2
# A[10,1:3] .= 1
# A[12,7] = 2
A[14,3] = 15
# A[14,3] = 5
# A[8,8] = 3
A[diagind(A)] .= collect(1:n)
# A[14,14] = 0
A[3,14] = 20

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

Y = zeros(n,nQueries)
for s in 1:nQueries
    Y[:,s] = A * X[:,s]
end

# p = heatmap2(X,Y)
# Plots.savefig(p, "/tmp/XY.png")

M = Convex.Variable(n,n)
# problem = minimize(norm(M,1), M*X == Y)
problem = minimize(sum(abs.(M)), M*X == Y)
# Convex.solve!(problem, () -> SCS.Optimizer(max_iters = 10000, verbose = false))
Convex.solve!(problem, () -> COSMO.Optimizer(verbose = false))
@show norm(evaluate(M),1)
@show norm(A,1)
# A*X ≈ evaluate(M)*X
# p = heatmap2(A,evaluate(M))
# Plots.savefig(p, "/tmp/XY.png")

pX = Plots.heatmap(X, title="X")
pY = Plots.heatmap(Y, title="Y")
cmin = minimum([A[:]..., evaluate(M)[:]...])
cmax = maximum([A[:]..., evaluate(M)[:]...])
pA = Plots.heatmap(A,clim=(cmin,cmax), title="A")
pAhat = Plots.heatmap(evaluate(M),clim=(cmin,cmax), title="Ahat")
cmin = minimum([(A*X)[:]..., (evaluate(M)*X)[:]...])
cmax = maximum([(A*X)[:]..., (evaluate(M)*X)[:]...])
pAX = Plots.heatmap(A*X,clim=(cmin,cmax), title="A*X")
pAhatX = Plots.heatmap(evaluate(M)*X,clim=(cmin,cmax), title="Ahat*X")
p = Plots.plot(pX, pA, pAX, pY,  pAhat, pAhatX, layout=(2,3))
Plots.savefig(p, "/tmp/not_optimal.png")
