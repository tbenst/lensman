using LinearAlgebra, StatsBase, ProgressMeter, JuMP, ANTsRegistration, Dates,
    NIfTI
import Base.Threads.@threads
sio = pyimport("scipy.io")

"""Sparse (min L1) solution to Y=AX for unknown A.
    
Pass ie `optimizer=Gurobi.Optimizer`
Recommend to use Gurobi, which is 50x faster than SCS.
Gurobi has a free academic license.
https://github.com/jump-dev/Gurobi.jl#installation
"""
function reconstructA(X,Y, optimizer)
    m,n = size(Y,1), size(X,1)
    model = Model(optimizer)
    @variable(model, Â[1:m,1:n])
    @variable(model, l1)
    @constraint(model, [l1; vec(Â)] in MOI.NormOneCone(m*n))
    @constraint(model, con, Â*X .== Y)
    @objective(model, Min, l1)
    optimize!(model)
    value.(Â)
end

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
    aa2a(stimGroups), aa2a(groupsPerCell)
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

"""
Add influence maps together for each stim that includes a particular neuron.
Return an entangled influence map for each neuron.

winSize: 3*volRate is reasonable for H2B-GCaMP6s
"""
function entangledInfluenceMaps(cells::DataFrame, tseries;
        winSize=10)
    roiMask = []
    (H, W, Z, T) = size(tseries)
    nCells = maximum(cells.cellIdx)
    entangledImaps = zeros(H, W, Z, nCells)

    @threads for c in 1:nCells
        stims = cells[cells.cellIdx .== c, [:stimStart, :stimStop]]
        iMap = zeros(H,W,Z)
        for (stimStart, stimStop) in eachrow(stims)
            f = mean(tseries[:,:,:,stimStop+1:stimStop+winSize], dims=4)[:,:,:,1]
            f0 = mean(tseries[:,:,:,stimStart-winSize:stimStart-1], dims=4)[:,:,:,1]
            df_f = @. (f - f0)/f0
            iMap .+= df_f
        end
        entangledImaps[:,:,:,c] .= iMap
    end
    
    entangledImaps
end

"Median filter, return same length vector."
function medianfilt(x, window=3)
    filtered = rollmedian(x, window)
    vcat([x[1]], filtered, [x[end]])
end

"""Use ANTs to perform registration, either affine or with SyN.

See https://academic.oup.com/view-large/120212390 for reference on using with
zebrafish atlas registration. That paper uses the following:

ants_register(fixed, moving; interpolation = "WelchWindowedSinc",
    histmatch = 0, sampling_frac = 0.25, maxiter = 200, threshold=1e-8,
    use_syn = false, synThreshold = 1e-7, synMaxIter = 200, dont_run = true)

Tyler previously used affine + sampling_frac=1.0 for faster online registration.

SyN took 9487s for a MultiMAP 

Writes intermediate data to tmp directory.
"""
function ants_register(fixed, moving; interpolation = "BSpline",
        histmatch = 1, maxiter = 200, threshold = "1e-6",
        initial_moving_type = 1, sampling_frac = 1.0,
        use_syn = false, synThreshold = 1e-7, synMaxIter = 200,
        ants_path = "/opt/ANTs/install/bin/antsRegistration",
        save_dir = ANTsRegistration.userpath(),
        dont_run = false)
    if use_syn
        mstr = "_SyN"
    else
        mstr = "_affine"
    end
    outprefix = Dates.format(DateTime(now()), DateFormat("YmmddTHHMMSSsZ"))
    fixedname = joinpath(save_dir, "$(outprefix)_fixed.nrrd")
    save(fixedname, fixed)
    movingname = joinpath(save_dir, "$(outprefix)_moving.nrrd")
    save(movingname, moving)
    outname = joinpath(save_dir, outprefix*mstr)
    # from ANTs zfish atlas paper
    if use_syn
        # syn_cmd = "-t SyN\[0.05,6,0.5\] -m CC\[$fixedname, $movingname,1,2\] -c \[$synMaxIter\x$synMaxIter\x$synMaxIter\x0,$synThreshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox"
        syn_cmd = `-t SyN\[0.05,6,0.5\] -m CC\[$fixedname, $movingname,1,2\] -c \[$synMaxIter\x$synMaxIter\x$synMaxIter\x$(synMaxIter)x10,$synThreshold,10\] --shrink-factors 12x8x4x2x1 --smoothing-sigmas 4x3x2x1x0vox`
    else
        syn_cmd = ``
    end
    affine_cmd = `$(ants_path) -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation $interpolation --use-histogram-matching $histmatch -r \[$fixedname, $movingname,$initial_moving_type\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,    $sampling_frac\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,$sampling_frac\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
    cmd = `$(affine_cmd) $(syn_cmd)`

    if dont_run
        return cmd
    end

    println("calling ANTs...")
    @time println(read(cmd, String))

    affine_transform_path = glob(outprefix*"*GenericAffine.mat", save_dir)[1]
    affine_transform = sio.loadmat(affine_transform_path)
    warpedname = joinpath(outname*"_Warped.nii.gz")

    # if this fails, tmpfs may be out of space
    # try `rm -rf /run/user/1000/tyler/ANTs/`
    registered = niread(warpedname);
    println("finished registering!")
    registered
end