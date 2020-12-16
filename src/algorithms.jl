using LinearAlgebra, StatsBase, ProgressMeter

num2String(n,base,nDigits) = join(reverse(digits(n-1, base=base, pad=nDigits)),"")

"Inverse of digits."
function undigits(d; base=10)
    (s, b) = promote(zero(eltype(d)), base)
    mult = one(s)
    for val in d
        s += val * mult
        mult *= b
    end
    return s
end

# we duplicate funcitons to support array of array (variable number cells
# per stim) or matrix (same number of cells per stim)
# array of array is 2x faster

function stonerStimGroups(nCells, base=2)
    nDigits = Int(ceil(log(base,nCells)))
    cells = []
    # represent each number as a vector ("string")
    # with proper base
    numVecs = [digits(c, base=base, pad=nDigits) for c in 0:nCells-1]
    stimGroups = []
    g = 0
    groupsPerCell = [[] for i in 1:nCells]
    for d in 1:nDigits
        for b in 0:base-1
            idxs = findall(getindex.(numVecs,d) .== b)
            g+=1
            push!(stimGroups, idxs)
            for c in idxs
                push!(groupsPerCell[c], g)
            end
        end
    end
    stimGroups, groupsPerCell
end

function stonerPerm(nCells, base=2)
    @assert base==2
    # otherwise, '111' -> '000' => duplicate string!!
    @assert log2(nCells) % 2 == 0
    nDigits = Int(ceil(log(base,nCells)))
    permutation = zeros(Int,nCells)
    for c in 1:nCells
        digitsArray = Bool.(digits(c-1,base=base, pad=nDigits))
        if sum(digitsArray) % 2 == 0
            # @show c, digitsArray
            permutation[c] = c
        else
            newDigitsArray = (~).(digitsArray)
            u = undigits(Int64.(newDigitsArray),base=2) + 1
            # @show c, u, digitsArray, newDigitsArray
            permutation[c] = u
        end
    end
    permutation
end

function repeatedStimGroups(nCells, base, nReps)
    
end


triangularNumber(n) = n*(n+1)/2

function l2OffDiagonal(mat)
    l2 = 0
    for i in 1:size(mat,1)
        for j in 1:size(mat,2)
            if i!=j
                l2 += mat[i,j]^2
            end
        end
    end
    l2
end

"""
Return concurrency score for StimGroups::Array{Array{Int,1}}

Each neuron i∈[1,N] is a member in  groups G = {g ∈ ℤᴺ | i ∈ g }


"""
function calc_concurrency_score(stimGroups)
    nCells = maximum(vcat(stimGroups...))
    score = 0
    for c in 1:nCells
        concurrencyCount = zeros(nCells)
        for group in stimGroups
            if c in group
                for d in group
                    concurrencyCount[d] += 1
                end
            end
        end
        deleteat!(concurrencyCount, c)
        # score += sum(triangularNumber.(concurrencyCount))
        score += norm(concurrencyCount)
    end
    score
end

function calc_concurrency_score(stimGroups::Array{T,2}) where T<:Real
    nCells = maximum(stimGroups)
    concurrencyCount = zeros(Int64, nCells, nCells)
    for c in 1:nCells
        for group in eachrow(stimGroups)
            if c in group
                for d in group
                    concurrencyCount[c,d] += 1
                end
            end
        end
    end
    l2OffDiagonal(concurrencyCount)
end

function count_concurrency(stimGroups)
    nCells = maximum(vcat(stimGroups...))
    counts = zeros(Int,nCells, nCells)
    for c in 1:nCells
        for group in stimGroups
            if c in group
                for d in group
                    counts[c,d] += 1
                end
            end
        end
    end
    counts
end

function count_concurrency(stimGroups::Array{T,2}) where T<:Real
    nCells = maximum(stimGroups)
    counts = zeros(Int,nCells, nCells)
    for c in 1:nCells
        for group in eachrow(stimGroups)
            if c in group
                for d in group
                    counts[c,d] += 1
                end
            end
        end
    end
    counts
end



function minMaxIntersection(stimGroups, groupsPerCell)
    # unique set intersection property
    maxIntersect = 1
    minIntersect = 1
    for groups in groupsPerCell
        setIntersection = intersect([stimGroups[g] for g in groups]...)
        maxIntersect = maximum([maxIntersect, length(setIntersection)])
        minIntersect = minimum([minIntersect, length(setIntersection)])
    end
    minIntersect, maxIntersect
end

function minMaxIntersection(stimGroups::Array{T,2}, groupsPerCell::Array{T,2}) where T<:Real
    # unique set intersection property
    maxIntersect = 1
    minIntersect = 1
    for groups in eachrow(groupsPerCell)
        setIntersection = intersect([stimGroups[g,:] for g in groups]...)
        maxIntersect = maximum([maxIntersect, length(setIntersection)])
        minIntersect = minimum([minIntersect, length(setIntersection)])
    end
    minIntersect, maxIntersect
end

function noDups(stimGroups)
    all(length.(stimGroups) .== length.(unique.(stimGroups)))
end

function noDups(stimGroups::Array{T,2}) where T<:Real
    all(length.(eachrow(stimGroups)) .== length.(unique.(eachrow(stimGroups))))
end


function calcGroupsPerCell(stimGroups, nCells, base)
    gpc = [[] for i in 1:nCells]
    for (g,group) in enumerate(stimGroups)
        for i in group
            push!(gpc[i], g)
        end
    end
    gpc
end

function calcGroupsPerCell(stimGroups::Array{T,2}, nCells, base) where T<:Real
    cellGroupCount = zeros(Int, nCells)
    for (g,group) in enumerate(eachrow(stimGroups))
        for c in group
            cellGroupCount[c] += 1
        end
    end
    # @show nGroupsPerCell = maximum(cellGroupCount)
    nGroupsPerCell = maximum(cellGroupCount)
    gpc = zeros(Int, nCells, nGroupsPerCell)
    groupIdx = ones(Int,nCells) # pointer to next idx to update per cell
    for (g,group) in enumerate(eachrow(stimGroups))
        for c in group
            gpc[c,groupIdx[c]] = g
            groupIdx[c] += 1
        end
    end
    gpc
end


function calcGroupsPerCell(stimGroups::Array{T,2}, nGroupsPerCell, nCells, base) where T<:Real
    gpc = zeros(Int, nCells, nGroupsPerCell)
    groupIdx = ones(Int,nCells) # pointer to next idx to update per cell
    for (g,group) in enumerate(eachrow(stimGroups))
        for c in group
            gpc[c,groupIdx[c]] = g
            groupIdx[c] += 1
        end
    end
    gpc
end

function printGroups(stimGroups,base,nDigits)
    println.(map(g->num2String.(g,base,nDigits), stimGroups))
end

function printGroups(stimGroups::Array{T,2},base,nDigits) where T<:Real
    println.(map(g->num2String.(g,base,nDigits), eachrow(stimGroups)))
end

function sameNstimsPerCell(groupsPerCell)
    all(length.(eachrow(groupsPerCell)) .== length(groupsPerCell[1,:]))
end

function randomSwaps(stimGroups, metric, nCells, base, nSwaps=10000;
  verbose=false, nDigits=4)
    stimGroups = deepcopy(stimGroups)
    proposedStimGroups = deepcopy(stimGroups)
    score = metric(stimGroups)
    @info "initial score: $score"
    @assert all(length.(stimGroups) .== length(stimGroups[1]))
    cellIdxs = collect(1:length(stimGroups[1]))
    groupIdxs = collect(1:length(stimGroups))
    if verbose
        println("Initial:\n")
        printGroups(stimGroups, base, nDigits)
    end
    for _ in 1:nSwaps
        a,b = sample(groupIdxs, (2,), replace=false)
        i,j = sample(cellIdxs, (2,))

        proposedStimGroups[a][i] = stimGroups[b][j]
        proposedStimGroups[b][j] = stimGroups[a][i]
        newScore = metric(proposedStimGroups)
        if (newScore < score)
            groupsPerCell = calcGroupsPerCell(proposedStimGroups, nCells, base)
            theMin, theMax = minMaxIntersection(proposedStimGroups, groupsPerCell)
            # TODO: test that num groups / cell are the same
            if theMin==1 & theMax==1 & noDups(proposedStimGroups)
                if verbose
                    println("Improved! New proposed:")
                    printGroups(proposedStimGroups, base, nDigits)
                end
                # return proposedStimGroups
                stimGroups = deepcopy(proposedStimGroups)
                score = newScore
            else
                proposedStimGroups = deepcopy(stimGroups)
            end
        else
             proposedStimGroups = deepcopy(stimGroups)
        end
    end
    @info "final score: $score"
    stimGroups
end


"More efficient when a 2D array (same length for all groups)."
function randomSwaps(stimGroups::Array{T,2}, groupsPerCell::Array{T,2}, metric, nCells, base, nSwaps=10000;
  verbose=false, nDigits=4) where T<:Real
    println("using efficient swap algo")
    stimGroups = deepcopy(stimGroups)
    proposedStimGroups = deepcopy(stimGroups)
    score = metric(stimGroups)
    @info "initial score: $score"
    @assert all(length.(eachrow(stimGroups)) .== length(stimGroups[1,:]))
    cellIdxs = collect(1:size(stimGroups, 2))
    groupIdxs = collect(1:size(stimGroups,1))
    
    nGroupsPerCell = size(groupsPerCell,2)
    
    if verbose
        println("Initial:\n")
        printGroups(stimGroups, base, nDigits)
    end
    @showprogress for _ in 1:nSwaps
        proposedStimGroups .= stimGroups
        # proposedStimGroups = deepcopy(stimGroups)
        a,b = sample(groupIdxs, (2,), replace=false)
        i,j = sample(cellIdxs, (2,))
        proposedStimGroups[a,i] = stimGroups[b,j]
        proposedStimGroups[b,j] = stimGroups[a,i]
        # println(noDups(proposedStimGroups), proposedStimGroups)
        newScore = metric(proposedStimGroups)
        if (newScore < score)
            # println("better score!", noDups(proposedStimGroups))
            groupsPerCell = calcGroupsPerCell(proposedStimGroups, nGroupsPerCell, nCells, base)
            theMin, theMax = minMaxIntersection(proposedStimGroups, groupsPerCell)
            # if theMin==1 & theMax==1 # many updates
            # if noDups(proposedStimGroups) #only one update?!
            if theMin==1 & theMax==1 & noDups(proposedStimGroups) & sameNstimsPerCell(groupsPerCell)
                if verbose
                    println("Improved! New proposed:")
                    printGroups(proposedStimGroups, base, nDigits)
                end
                # return proposedStimGroups
                stimGroups .= proposedStimGroups
                score = newScore
            end
        end
    end
    @info "final score: $score"
    groupsPerCell = calcGroupsPerCell(stimGroups, nGroupsPerCell, nCells, base)
    stimGroups, groupsPerCell
end