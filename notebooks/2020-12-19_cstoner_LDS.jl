using Lensman, Random, Gadfly, DataFrames, StatsBase
import Cairo, Fontconfig
using Plots, LinearAlgebra, PyCall
using Convex, SCS

py"""
import numpy as np
import scipy

def makeStableMatrix(n, density=0.1):
    # TODO: implement on GPU?
    eig = np.diag(np.random.uniform(-1,0,n))
    v = scipy.linalg.orth(np.random.uniform(-1,1,(n,n)))
    A = v.T @ eig @ v
    # need to make A sparse
    indices = np.reshape(np.indices((n,n)),[2,-1])
    mask = np.eye(n)
    # may create density slightly less than desired if choice is on diagonal
    idx = indices[:,np.random.choice(np.arange(indices.shape[1]),int(np.ceil(density*n**2-n)),replace=False)]
    mask[idx[0], idx[1]] = 1
    A[np.logical_not(mask)] = 0
    return A
"""
makeStableMatrix(n,density=0.1) = py"makeStableMatrix"(n, density)
##

N = 256
base = 2
stimGroups, groupsPerCell = stonerStimGroups(N, base)
stimGroupsArray = copy(hcat(stimGroups...)')
groupsPerCellArray = copy(hcat(groupsPerCell...)')

stoner_perm = stonerPerm(N)
stimGroupsPerm = map(g->perm(stoner_perm,g),stimGroups)
stonerGroups = aa2a(vcat(stimGroups, stimGroupsPerm))
stonerGroupsPerCell = calcGroupsPerCell(stonerGroups, N, base)


optGroups, optGpc = randomSwaps(stimGroupsArray, groupsPerCellArray, calc_concurrency_score,
N, base, 100)

optStonerGroups, optStonerGpc = randomSwaps(stonerGroups, stonerGroupsPerCell, calc_concurrency_score,
N, base, 100)

concurrencyMat = count_concurrency(stimGroups)
optConcurrencyMat = count_concurrency(optGroups)
stonerMat = count_concurrency(stonerGroups)
optStonerMat = count_concurrency(optStonerGroups)

C = concurrencyMat
stimGroups = stimGroupsArray
groupsPerCell = groupsPerCellArray

# C = stonerMat
# stimGroups = stonerGroups
# groupsPerCell = stonerGroupsPerCell

# C = optConcurrencyMat
# stimGroups = optGroups
# groupsPerCell = optGpc

C = optStonerMat
stimGroups = optStonerGroups
groupsPerCell = optStonerGpc

@show rank.([concurrencyMat, stonerMat, optConcurrencyMat, optStonerMat])
# A = makeStableMatrix(N,1.)
A = zeros(N,N)
# A = reshape(shuffle(collect(1:N*N)),N,N)
# A[sample(CartesianIndices(A)[:], 1, replace=false)] .= 1
A[10,50] = 1
heatmap(A)

# OLS
T = 1000*N
Xall = zeros(N,T)
Xall[:,1] .= randn(N)
sigma = 0.001
for t in 2:T
    Xall[:,t] .= A*Xall[:,t-1] + sigma .* randn(N)
end
X = Xall[:,1:end-1]
Y = Xall[:,2:end]

A_ols = Y * transpose(X) * inv(X * transpose(X))
# A_ols = transpose(X) \ transpose(Y)
# A_ols = transpose(Y) \ transpose(X)
#
# p_ols = heatmap(A_ols)
# p_A = heatmap(A)
# p = heatmap(A-A_ols)
# l = @layout [a b]
# Plots.plot(p_A,p_ols, layout=l)


# combo stim
nStim = size(stimGroups,1)

XcAll = zeros(N,nStim*2)
stimMagnitude = 1
for s in 1:nStim
    neuronsToStim = stimGroups[s,:]
    XcAll[neuronsToStim,s*2-1] .= stimMagnitude
end
for t in 2:2:nStim*2
    XcAll[:,t] = A * XcAll[:,t-1]
end

entangledInfluenceMaps = []
for (c,groups) in enumerate(eachrow(groupsPerCell))
    entangledInfluenceMap = sum(XcAll[:,2 .* groups],dims=2)
    push!(entangledInfluenceMaps, entangledInfluenceMap)
end

entangledInfluenceArray = hcat(entangledInfluenceMaps...)
Gadfly.plot(y=sum(entangledInfluenceArray,dims=1))
# as expected, only #50 is 0...
# so why don't we reconstruct...?
Gadfly.plot(y=[sum(abs.(x .- entangledInfluenceArray[10,:])) for x in eachrow(concurrencyMat)])
##
# Î = pinv(C) * transpose(entangledInfluenceArray)
# add l1..?
# Î = pinv(C) * transpose(entangledInfluenceArray)

# CI = fancyI
all(A*C .== entangledInfluenceArray)
all(C*transpose(A) .== transpose(entangledInfluenceArray))
all(transpose(A) .== pinv(C) * transpose(entangledInfluenceArray))

p_A = heatmap(A,title="A")
p_entang = heatmap(entangledInfluenceArray,title="entangled")
p_C = heatmap(C,title="concurrencyMat")
p_Ihat = heatmap(transpose(evaluate(Î)),title="A reconstructed")
p = Plots.plot(p_A,p_entang, p_C,p_Ihat, layout=4)
# savefig(p,"/tmp/single_entry.png")

# heatmap(C*Î)

# heatmap(entangledInfluenceArray .- concurrencyMat)

## minimum / percentile...

minInfluenceMaps = []
for (c,groups) in enumerate(eachrow(groupsPerCell))
    minInfluenceMap = minimum(XcAll[:,2 .* groups],dims=2)
    push!(minInfluenceMaps, minInfluenceMap)
end
minInfluenceMaps = copy(transpose(aa2a(minInfluenceMaps)))
p_A = heatmap(A,title="A")
p_mI = heatmap(minInfluenceMaps)
p = Plots.plot(p_A,p_mI, layout=2)
