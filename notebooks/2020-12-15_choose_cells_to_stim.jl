##
ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random
import Gadfly
using Unitful: μm, m, s
##
# tif840path = "/mnt/deissero/users/tyler/b115/2020-12-15_h2b6s_chrmine_kv2.1_5dpf/fish2/SingleImage-840nm-512-016/SingleImage-840nm-512-016_Cycle00001_Ch3_000001.ome.tif"
tif840path = "/mnt/deissero/users/tyler/b115/2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2/SingleImage-840nm-1024-017/SingleImage-840nm-1024-017_Cycle00001_Ch3_000001.ome.tif"
## to burn at etl=0 if using calibration circa 2020-12-15, need +45 offset
offset = float(uconvert(m, 45μm)) / m
fishDir = joinpath(splitpath(tif840path)[1:end-2]...)
expName = splitpath(tif840path)[end-1]
avgImage = ImageMagick.load(tif840path)
size(avgImage)
##
# compile
##
candNeuronCenterMask = findNeurons(avgImage,1.5,2,10);
candidateTargetLocs = findall(candNeuronCenterMask)
# candidateTargets = copy(mask2ind(candNeuronCenterMask))

# addTargetsToImage(avgImage, candidateTargetLocs)


img = RGB.(imadjustintensity(avgImage))
stim_points = zeros(Bool,size(img))
stim_points[candidateTargetLocs] .= true
stim_points = dilate(dilate(stim_points))
channelview(img)[[1,3],:,:,:] .= 0
channelview(img)[1,:,:,:] .= 0.25 * float(stim_points)
imshow(img)

## Sample 128 neurons
nNeurons = 1024
# nNeurons = 1
neuron_locs = sample(candidateTargetLocs, nNeurons, replace=false)     

# Visualize stim targets
im = RGB.(imadjustintensity(avgImage))
stim_points = zeros(Bool,size(im))
stim_points[neuron_locs] .= true
stim_points = dilate(dilate(stim_points))
channelview(im)[[1,3],:,:,:] .= 0
channelview(im)[1,:,:,:] .= 0.5*float(stim_points)
# imshow(im)
im;

##
# targets1 = [ 606. 636. offset;
#            ] ./ 2

# create_slm_stim([targets1],
#     "/mnt/deissero/users/tyler/slm/masks/2020-11-30/fish1-single-cell")


# one cell at a time
# target_groups = [cartIdx2SeanTarget(neuron_loc, offset)
#     for neuron_loc in neuron_locs]

# k neurons per stim
k = 32
target_groups = [vcat(cartIdx2SeanTarget.(locs, fill(offset, k))...)
    for locs in Iterators.partition(neuron_locs,k)]


## Save files for SLM stim
name = "1024cell-32concurrent-zoffset_fix"
outname = joinpath(fishDir, name)
save("$(outname).tif", im)

create_slm_stim(target_groups,
    joinpath(fishDir, name))



## Then use image-block ripping, and we open 
# TODO replace with streaming...?


############# Analysis ##################
## load in results
# only triggered the first 120 trials, oops!!
fishDir = "/media/tyler/tyler_benster/2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2"
tseriesDir = joinpath(fishDir, "TSeries-1024cell-32concurrent-zoffset_fix-043")
slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"

tseries = loadTseries(tseriesDir);
(H, W, Z, T) = size(tseries)
##
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
allWithin(diff(stimStartIdx),0.05)
# adjust for expected number of stimuli....
@assert length(stimStartIdx) == 32
allWithin(stimEndIdx .- stimStartIdx,0.05)
avgStimDuration = mean(stimEndIdx .- stimStartIdx)
# adjust for volume rate / stim duration...
@assert 30 > avgStimDuration > 3
# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

@warn "may not be needed if in same session / will clobber..."
# Assumes no sequence of stim
target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
##
@assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}

targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
    is1024=false, zOffset=45)

cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)

stimLocs = map(CartesianIndex ∘ Tuple, eachrow(cells[:,[2,1]]))
first(cells)
##
stimPlane = Int(mean(map(x->mean(x[:,3]), targetsWithPlaneIndex)))
avgImage = dropdims(mean(tseries[:,:,stimPlane,1:1000], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
avgImageAdj = RGB.(avgImageAdj)
channelview(avgImageAdj)[[1,3],:,:,:] .= 0
@warn "hardcode for 512x512"
microscope_units = (2* 0.6299544139175637μm, 2* 0.6299544139175637μm, 2.0μm)
# assume that "7um" spiral galvo signal is calibrated for 16x, not 25x
# "16x" is actually 14.4x
targetSizePx = 7μm * (14.4/25) / microscope_units[1]
addTargetsToImage(copy(avgImageAdj), cartIdx2Array(stimLocs),
    targetSize=targetSizePx)
##
cellDf_f = []
# @assert size(stimStartIdx,1)==size(stimLocs,1) # only did 120 by accident
winSize = 15
# TODO: don't hardcode frameRate (15)
# delay=Int(ceil(15*2))
delay=0
roiMask = []
# TODO: assumes one z plane!
for (x,y,z,stimStart,stimStop) = eachrow(cells)
    mask = Gray.(zeros(Bool, H,W,Z))
    draw!(view(mask,:,:,z), Ellipse(CirclePointRadius(x,y,targetSizePx)))
    mask = findall(channelview(mask))
    push!(roiMask, mask)
    # hardcode 1 plane 
    neuronTrace = extractTrace(tseries, mask)
    neuronTrace = imageJkalmanFilter(neuronTrace)
    f = mean(neuronTrace[stimStop+delay:stimStop+winSize])
    
    # problem: if neuron happens to be active and declining before stim,
    # we may miss
    f0 = mean(neuronTrace[stimStart-winSize:stimStart-1])

    # solution: use percentile f0
    # f0 = percentile(neuronTrace,0.1)
    df_f = (f-f0)/(f0+1e-8)
    # df_f = (f-f0)
    push!(cellDf_f,df_f)
end

if "df_f" in names(cells)
    cells.df_f = cellDf_f
else
    insertcols!(cells, size(cells,2)+1, :df_f => cellDf_f)
end

rankings = sortperm(cellDf_f)
best = rankings[end]
worst = rankings[1]

idx = rankings[end-10]

sum(cellDf_f.>0.4)
# imshow(addTargetsToImage(avgImage, stimLocs[[best],:]))
# imshow(tseries[:,:,1,stimStartIdx[idx]-winSize:stimEndIdx[idx]+winSize])
##
cmax = maximum(abs.(cellDf_f))
p_df_map = Gadfly.with_theme(:dark) do
    Gadfly.plot(cells, x=:x, y=:y, Gadfly.Geom.point, color=:df_f,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true))
end

img = PNG(joinpath(plotDir, "$(expName)_df_f_map.png"), 6inch, 4inch)
p_df_map
## viz top 64 (blue) vs worst 64 (red)
img = RGB.(avgImageAdj)
channelview(img)[[1,3],:,:,:] .= 0
img = addTargetsToImage(img, cartIdx2Array(stimLocs[rankings[end-63:end]]),
    channel="blue", targetSize=targetSizePx)
img = addTargetsToImage(img, cartIdx2Array(stimLocs[rankings[1:64]]),
    channel="red",alpha=0.6, targetSize=targetSizePx)
img


##
function plotStim(tseries,roiMask,cells, idx)
    roiM = roiMask[idx]
    x, y, z, stimStart, stimStop = cells[idx,:]
    plotRange = stimStart-30:stimStop+30
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    fluorescentTrace = imageJkalmanFilter(fluorescentTrace)
    p = plot(plotRange/frameRate, fluorescentTrace,
        xlabel="time (s)", ylabel="raw fluorescence")
    plot!(Lensman.rectangle((stimStop-stimStart)/frameRate,maximum(fluorescentTrace),
        stimStart/frameRate,0), opacity=.5, label="stim")
    p
end
# best
plots = []
nplots = 64
for idx in rankings[end:-1:end-nplots+1]
    push!(plots, plotStim(tseries,roiMask,cells,idx))
end
p = plot(plots..., layout=(16, 4), legend=false, size=(870,1024*2))
plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

savefig(p, joinpath(plotDir, "$(expName)_best64traces.png"))
p
##
plots = []
nplots = 64
plotIdxOffset = 575
for idx in rankings[end-plotIdxOffset:-1:end-plotIdxOffset-nplots+1]
    push!(plots, plotStim(tseries,roiMask,cells,idx))
end
plot(plots..., layout=(16, 4), legend=false, size=(870,1024*2))

##
Gadfly.plot(cells, x=:df_f, Gadfly.Geom.histogram)


### GENERATE COMBO Experiment


# select cells with df/f >= 1.0
min_df_F = 1.0
maxIdx = searchsortedfirst(sort(cells, :df_f, rev=true).df_f, min_df_F, rev=true)
selectedStimLocs = sample(stimLocs[rankings[end:-1:end-maxIdx+2]], 64, replace=false)

##
nCells = 64
base = 2
nReps = 8
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


stonerGroups
##
k = Int(nCells / base)
target_groups = []
@warn "512"
for group in eachrow(stonerGroups)
    groupStimLocs = selectedStimLocs[group]
    push!(target_groups, vcat(cartIdx2SeanTarget512.(groupStimLocs, fill(offset, k))...))
end

## Save files for SLM stim
slmMaskOutputDir = joinpath(splitpath(tif840path)[1:end-2]...)

name = "cstoner_n64_b2_r8"
outname = joinpath(slmMaskOutputDir, name)
# save("$(outname).tif", im)

create_slm_stim(target_groups,
    joinpath(slmMaskOutputDir, name))

##
addTargetsToImage(copy(avgImageAdj), cartIdx2Array(selectedStimLocs),
    targetSize=targetSizePx)
