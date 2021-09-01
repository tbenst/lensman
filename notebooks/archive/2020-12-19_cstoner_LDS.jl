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

N = 100
base = 10
stimGroups, groupsPerCell = stonerStimGroups(N, base)
stimGroupsArray = copy(hcat(stimGroups...)')
groupsPerCellArray = copy(hcat(groupsPerCell...)')

# stoner_perm = stonerPerm(N)
# stimGroupsPerm = map(g->perm(stoner_perm,g),stimGroups)
# stonerGroups = aa2a(vcat(stimGroups, stimGroupsPerm))
# stonerGroupsPerCell = calcGroupsPerCell(stonerGroups, N, base)


# optGroups, optGpc = randomSwaps(stimGroupsArray, groupsPerCellArray, calc_concurrency_score,
# N, base, 100)

# optStonerGroups, optStonerGpc = randomSwaps(stonerGroups, stonerGroupsPerCell, calc_concurrency_score,
# N, base, 100)

concurrencyMat = count_concurrency(stimGroups)
# optConcurrencyMat = count_concurrency(optGroups)
# stonerMat = count_concurrency(stonerGroups)
# optStonerMat = count_concurrency(optStonerGroups)

# @show rank.([concurrencyMat, stonerMat, optConcurrencyMat, optStonerMat])

##

C = concurrencyMat
stimGroups = stimGroupsArray
groupsPerCell = groupsPerCellArray

# C = stonerMat
# stimGroups = stonerGroups
# groupsPerCell = stonerGroupsPerCell

# C = optConcurrencyMat
# stimGroups = optGroups
# groupsPerCell = optGpc

# C = optStonerMat
# stimGroups = optStonerGroups
# groupsPerCell = optStonerGpc
##
# A = makeStableMatrix(N,0.1)
A = reshape(shuffle(collect(1:N*N)),N,N)
A[StatsBase.sample(CartesianIndices(A)[:], N*N - N, replace=false)] .= 0
# A = zeros(N,N)
# A[10,50] = 1
# heatmap(A)
##
# OLS
# T = 1000*N
# Xall = zeros(N,T)
# Xall[:,1] .= randn(N)
# sigma = 0.001
# for t in 2:T
#     Xall[:,t] .= A*Xall[:,t-1] + sigma .* randn(N)
# end
# X = Xall[:,1:end-1]
# Y = Xall[:,2:end]

# A_ols = Y * transpose(X) * inv(X * transpose(X))
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

# the simple intersection-of-group strategy
minInfluenceMaps = []
for (c,groups) in enumerate(eachrow(groupsPerCell))
    # Xgroup = XcAll[:,2 .* groups]
    # minIdxs = argmin(abs.(Xgroup), dims=2)
    # minInfluenceMap = sign.(Xgroup[minIdxs]) .* Xgroup[minIdxs]
    minInfluenceMap = minimum(XcAll[:,2 .* groups],dims=2)
    push!(minInfluenceMaps, minInfluenceMap)
end
minInfluenceMaps = copy(transpose(aa2a(minInfluenceMaps)))

# # the entangled strategy
entangledInfluenceMaps = []
for (c,groups) in enumerate(eachrow(groupsPerCell))
    entangledInfluenceMap = sum(XcAll[:,2 .* groups],dims=2)
    push!(entangledInfluenceMaps, entangledInfluenceMap)
end
entangledInfluenceArray = hcat(entangledInfluenceMaps...)

# Gadfly.plot(y=sum(entangledInfluenceArray,dims=1))
# as expected, only #50 is 0...
# so why don't we reconstruct...?
# Gadfly.plot(y=[sum(abs.(x .- entangledInfluenceArray[10,:])) for x in eachrow(concurrencyMat)])

# Î = pinv(C) * transpose(entangledInfluenceArray)
# add l1..?
# Î = pinv(C) * transpose(entangledInfluenceArray)

# CI = fancyI
# all(A*C .== entangledInfluenceArray)
# all(C*transpose(A) .== transpose(entangledInfluenceArray))
# all(transpose(A) .== pinv(C) * transpose(entangledInfluenceArray))
A0Idx = A .== 0
# println("frac of matching: $(sum(A .== minInfluenceMaps)/length(A))")
# println("frac of zero matching: $(sum(A[A0Idx] .== minInfluenceMaps[A0Idx])/sum(A0Idx))")
# println("frac of non-zero matching: $(sum(A[(~).(A0Idx)] .== minInfluenceMaps[(~).(A0Idx)])/sum((~).(A0Idx)))")

M1Idx = minInfluenceMaps .!= 0
# println("frac of matching: $(sum(A .== minInfluenceMaps)/length(A))")
# println("Type II error rate: $(1 - sum(A[A0Idx] .== minInfluenceMaps[A0Idx])/sum(A0Idx))")
# println("Type I error rate: $(1 - sum(A[M1Idx] .== minInfluenceMaps[M1Idx])/sum(M1Idx))")

true_positives = sum(A[(~).(A0Idx)] .== minInfluenceMaps[(~).(A0Idx)])
false_negatives = sum(A[(~).(A0Idx)] .!= minInfluenceMaps[(~).(A0Idx)])
true_negatives = sum(A[A0Idx] .== minInfluenceMaps[A0Idx])
false_positives = sum(A[A0Idx] .!= minInfluenceMaps[A0Idx])

println("Precision: $(true_positives/(true_positives+false_positives))")
println("recall: $(true_positives/(true_positives+false_negatives))")
println("Type I error rate: $(1 - sum(A[M1Idx] .== minInfluenceMaps[M1Idx])/sum(M1Idx))")

##
p_A = heatmap(A,title="A")
p_mI = heatmap(minInfluenceMaps)
p = Plots.plot(p_A,p_mI, layout=2)

# p_A = heatmap(A,title="A")
# p_entang = heatmap(entangledInfluenceArray,title="entangled")
# p_C = heatmap(C,title="concurrencyMat")
# p_Ihat = heatmap(transpose(evaluate(Î)),title="A reconstructed")
# p = Plots.plot(p_A,p_entang, p_C,p_Ihat, layout=4)
# savefig(p,"/tmp/single_entry.png")

# heatmap(C*Î)

# heatmap(entangledInfluenceArray .- concurrencyMat)

## minimum / percentile...


p_A = heatmap(A,title="A")
p_mI = heatmap(minInfluenceMaps)
p = Plots.plot(p_A,p_mI, layout=2)

## Lasso
using ScikitLearn, COSMO
@sk_import linear_model: Lasso

lm = ScikitLearn.fit!(Lasso(alpha=0.1), C, entangledInfluenceArray)
println(lm.coef_)

## convex playground
X = Convex.Variable(N, N)
B = entangledInfluenceArray
A = C
λ = 10
problem = minimize(norm(A*X - B,2),[norm(X,1) <= λ])
Convex.solve!(problem, () -> COSMO.Optimizer(verbose = false))
