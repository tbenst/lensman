using Lensman, Random, Gadfly, DataFrames
import Cairo, Fontconfig
# debug / warmup

perm(permutation, array) = map(x->permutation[x],array)
nCells = 8
base = 2
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
stimGroupsArray = copy(hcat(stimGroups...)')
groupsPerCellArray = copy(hcat(groupsPerCell...)')

optGroups, optGpc = randomSwaps(stimGroupsArray, groupsPerCellArray, calc_concurrency_score,
    nCells, base, 100,verbose=true)


count_concurrency(optGroups)
count_concurrency(stimGroups)
## benchmark
nCells = 8
base = 2
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
@btime optGroups = randomSwaps(stimGroups, calc_concurrency_score,
    nCells, base, 100000)

# 6x speedup! :)
stimGroupsArray = copy(hcat(stimGroups...)')
groupsPerCellArray = copy(hcat(groupsPerCell...)')
@btime optGroups = randomSwaps(stimGroupsArray, groupsPerCellArray, calc_concurrency_score,
    nCells, base, 100000)
##


nCells = 64
base = 2
nReps = 8
println("==== ($base, $nCells, $nReps) ====")
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
stimGroups, groupsPerCell = aa2a.([stimGroups, groupsPerCell])
sameGroups = stimGroups
permGroups = stimGroups
stonerGroups = stimGroups
for i in 1:nReps-1
    permutation = randperm(nCells)
    if i % 2 == 1
        stoner_perm = stonerPerm(nCells)
    else
        stoner_perm = collect(1:nCells)
    end
    stimGroupsPerm = map(g->perm(permutation,g),stimGroups)
    stimGroupsPerm2 = map(g->perm(stoner_perm,g),stimGroups)
    sameGroups = vcat(sameGroups, stimGroups)
    permGroups = vcat(permGroups, stimGroupsPerm)
    stonerGroups = vcat(stonerGroups, stimGroupsPerm2)
end

sameGroupsPerCell = calcGroupsPerCell(sameGroups, nCells, base)
permGroupsPerCell = calcGroupsPerCell(permGroups, nCells, base)
stonerGroupsPerCell = calcGroupsPerCell(stonerGroups, nCells, base)

@assert count_concurrency(sameGroups)[1,1] == size(groupsPerCell,2)*nReps
@assert sum(count_concurrency(sameGroups)) == sum(count_concurrency(permGroups))
@assert sum(count_concurrency(sameGroups)) == sum(count_concurrency(stonerGroups))
# sameGroupsOpt1, _ =  randomSwaps(sameGroups, sameGroupsPerCell, calc_concurrency_score,
#     nCells, base, 1000)

permGroupsOpt, permGroupsPerCellOpt = randomSwaps(permGroups, permGroupsPerCell, calc_concurrency_score,
nCells, base, 10000)
# count_concurrency(permGroupsPerCellOpt)

stonerGroupsOpt, stonerGroupsPerCellOpt = randomSwaps(stonerGroups, stonerGroupsPerCell, calc_concurrency_score,
nCells, base, 10000)


@assert sum(count_concurrency(sameGroups)) == sum(count_concurrency(permGroupsOpt))
@assert sum(count_concurrency(sameGroups)) == sum(count_concurrency(stonerGroupsOpt))
@assert count_concurrency(sameGroups)[1,1] == count_concurrency(stonerGroupsOpt)[1,1]
@show calc_concurrency_score(sameGroups)
@show calc_concurrency_score(permGroups)
@show calc_concurrency_score(stonerGroups)
@show calc_concurrency_score(permGroupsOpt)
@show calc_concurrency_score(stonerGroupsOpt)
# @show calc_concurrency_score(sameGroupsOpt1)
# @show calc_concurrency_score(permGroupsOpt2)
# @show calc_concurrency_score(stonerGroupsOpt)

# count_concurrency(stimGroups)
# count_concurrency(stonerGroups)
# count_concurrency(stimGroupsOpt)
# count_concurrency(permGroupsOpt2)
# count_concurrency(stonerGroupsOpt)

# stoner even/odd gives 960 x #8 in ncells=64,base=2,nReps=2
# after 10000 swaps, we have only 230 x #8 BUT we have 6 x #9 !!

c_same = count_concurrency(sameGroups)
@show c_same_max = maximum(c_same[c_same.<c_same[1,1]])
sum(c_same.==c_same_max)
sum(c_same.==c_same_max - 1)

c_stoner = count_concurrency(stonerGroups)
@show c_stoner_max = maximum(c_stoner[c_stoner.<c_stoner[1,1]])
sum(c_stoner.==c_stoner_max)
sum(c_stoner.==c_stoner_max - 1)

c_perm = count_concurrency(permGroups)
@show c_perm_max = maximum(c_perm[c_perm.<c_perm[1,1]])
sum(c_perm.==c_perm_max)
sum(c_perm.==c_perm_max - 1)

c_permOpt = count_concurrency(permGroupsOpt)
@show c_permOpt_max = maximum(c_permOpt[c_permOpt.<c_permOpt[1,1]])
sum(c_permOpt.==c_permOpt_max)
sum(c_permOpt.==c_permOpt_max - 1)

c_stonerOpt = count_concurrency(stonerGroupsOpt)
@show c_stonerOpt_max = maximum(c_stonerOpt[c_stonerOpt.<c_stonerOpt[1,1]])
sum(c_stonerOpt.==c_stonerOpt_max)
sum(c_stonerOpt.==c_stonerOpt_max - 1)

# # better!!
# c_opt = count_concurrency(permGroupsOpt2)
# @show c_opt_max = maximum(c_opt[c_opt.<c_opt[1,1]])
# sum(c_opt.==c_opt_max)
# sum(c_opt.==c_opt_max-1)

# c_stoner_opt = count_concurrency(stonerGroupsOpt)
# @show c_stoner_opt_max = maximum(c_stoner_opt[c_stoner_opt.<c_stoner_opt[1,1]])
# sum(c_stoner_opt.==c_stoner_opt_max)
# sum(c_stoner_opt.==c_stoner_opt_max-1)


# count_concurrency(permGroupsPerCell)

## randperm on c_stoner only


# nCells = 256
# base = 2
# nReps = 8
nCells = 64
# nCells = 16
base = 2
nReps = 2
# nCells = 64
# base = 2
# nReps = 8
println("==== ($base, $nCells, $nReps) ====")
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
# count_concurrency(stimGroups)
#
stimGroups, groupsPerCell = aa2a.([stimGroups, groupsPerCell])
stoner_perm = stonerPerm(nCells)
stonerGroupsPerm = map(g->perm(stoner_perm,g),stimGroups)
stonerGroupsBase = vcat(stimGroups, stonerGroupsPerm)
stonerGroups = copy(stonerGroupsBase)

for i in 2:Int(nReps/2)
    permutation = randperm(nCells)
    tempStimGroupsPerm = map(g->perm(permutation,g),stonerGroupsBase)
    stonerGroups = vcat(stonerGroups, tempStimGroupsPerm)
    # stonerGroups = vcat(stonerGroups, stonerGroupsBase)
end

stonerGroupsPerCell = calcGroupsPerCell(stonerGroups, nCells, base)
stonerGroupsOpt, stonerGroupsPerCellOpt = randomSwaps(stonerGroups, stonerGroupsPerCell, calc_concurrency_score,
nCells, base, 10000)


c_stoner = count_concurrency(stonerGroups)
@show c_stoner_max = maximum(c_stoner[c_stoner.<c_stoner[1,1]])
@show sum(c_stoner.==c_stoner_max)
@show sum(c_stoner.==c_stoner_max - 1)
c_stoner


c_stonerOpt = count_concurrency(stonerGroupsOpt)
@show c_stonerOpt_max = maximum(c_stonerOpt[c_stonerOpt.<c_stonerOpt[1,1]])
@show sum(c_stonerOpt.==c_stonerOpt_max)
@show sum(c_stonerOpt.==c_stonerOpt_max - 1)

c_stonerOpt
df = DataFrame(concurrency_count=c_stoner[:], method="c_stoner")
df = vcat(df, DataFrame(concurrency_count=c_stonerOpt[:], method="c_stonerOpt"))
nbins = maximum(c_stoner[:])+1

grid = Geom.subplot_grid(Geom.histogram())
p = plot(df, x=:concurrency_count, grid, ygroup=:method, Scale.y_log2())
# draw(PNG("/home/tyler/Downloads/gadfly_hist_log2.png", 4inch, 3inch), p)
