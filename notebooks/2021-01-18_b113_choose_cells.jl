##
ENV["DISPLAY"] = "localhost:10.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, LightXML
import Gadfly
using Unitful: μm, m, s
##
tseriesDir = "/mnt/deissero/users/tyler/b113/2021-01-18_elavl3-chrmine-kv2.1_h2b6s_6dpf/fish1_chrmine/resting-005"
tseries = loadTseries(tseriesDir)

fishDir = joinpath(splitpath(tseriesDir)[1:end-2]...)
expName = splitpath(tseriesDir)[end-1]
avgImage = mean(tseries, dims=(3,4))[:,:,1,1]
size(avgImage)
imshow(avgImage)
## may want to now skip to Analysis section if not analyzing live
# compile
##
candNeuronCenterMask = findNeurons(avgImage,thresh_adjust=1.5, featSize=4,maxiSize=8);
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
println("found $(size(candidateTargetLocs,1)) potential targets")
## Sample 128 neurons
nNeurons = 1024
# not needed if not so many...
if size(candidateTargetLocs,1) < nNeurons
    neuron_locs = candidateTargetLocs
    nNeurons = size(candidateTargetLocs,1)
else
    neuron_locs = sample(candidateTargetLocs, nNeurons, replace=false)
end

# Visualize stim targets
img = RGB.(imadjustintensity(avgImage))
stim_points = zeros(Bool,size(img))
stim_points[neuron_locs] .= true
stim_points = dilate(dilate(stim_points))
channelview(img)[[1,3],:,:,:] .= 0
channelview(img)[1,:,:,:] .= 0.5*float(stim_points)
imshow(img)
img;

## Write xml
# create an empty XML document
W = size(avgImage,2) 
xmlGpl = XMLDocument() # create & attach a root node
xGplRoot = create_root(xmlGpl, "PVGalvoPointList")
xmlSeries = XMLDocument() # create & attach a root node
xSeriesRoot = create_root(xmlSeries, "PVSavedMarkPointSeriesElements")
set_attribute(xSeriesRoot, "Iterations", "1")
set_attribute(xSeriesRoot, "IterationDelay", "5000.0")
for (i,ci) in enumerate(neuron_locs[randperm(nNeurons)])
    y,x = Tuple(ci)
    # empirically, X ∈ [-7.6, 7.6], Y ∈ [8.3, -8.3]
    # where "(1,1)" is (-7.6,8.3)
    # and "(512,512)" is (7.6,-8.3)
    # and "(1,74)" meaning row 1 col 74 is (-5.37,8.3) 
    # and "(256,256)" is (0,0)
    x = (x-W/2)/(W/2) # map to (-1,1)
    x *= 7.6
    y = (y-W/2)/(W/2)
    y *= 8.3
    
    gplList = new_child(xGplRoot, "PVGalvoPoint")
    # todo use Dict or named argument instead..?
    set_attribute(gplList, "X", "$x")
    set_attribute(gplList, "Y", "$y")
    set_attribute(gplList, "Name", "Point $i")
    set_attribute(gplList, "Index", "$(i-1)")
    set_attribute(gplList, "ActivityType", "MarkPoints")
    set_attribute(gplList, "UncagingLaser", "Uncaging")
    set_attribute(gplList, "UncagingLaserPower", "0.76")
    set_attribute(gplList, "Duration", "2")
    set_attribute(gplList, "IsSpiral", "True")
    set_attribute(gplList, "SpiralSize", "0.220801205703228") # 7μm
    set_attribute(gplList, "SpiralRevolutions", "5")
    set_attribute(gplList, "Z", "0")

    mpe = new_child(xSeriesRoot, "PVMarkPointElement")
    set_attribute(mpe, "Repetitions", "5")
    set_attribute(mpe, "UncagingLaser", "Uncaging")
    set_attribute(mpe, "UncagingLaserPower", "0")
    set_attribute(mpe, "TriggerFrequency", "None")
    set_attribute(mpe, "TriggerSelection", "None")
    set_attribute(mpe, "TriggerCount", "1")
    set_attribute(mpe, "AsyncSyncFrequency", "None")
    set_attribute(mpe, "VoltageOutputCategoryName", "None")
    set_attribute(mpe, "VoltageRecCategoryName", "Current")
    set_attribute(mpe, "parameterSet", "CurrentSettings")

    # nested under <PVMarkPointElement>!
    gpe = new_child(mpe, "PVGalvoPointElement")
    if i==1
        set_attribute(gpe, "InitialDelay", "15")
    else
        set_attribute(gpe, "InitialDelay", "0.12")
    end
    set_attribute(gpe, "InterPointDelay", "248")
    set_attribute(gpe, "Duration", "170") # from exported file but ????
    # set_attribute(gpe, "Duration", "2")
    set_attribute(gpe, "SpiralRevolutions", "5")
    set_attribute(gpe, "AllPointsAtOnce", "False")
    set_attribute(gpe, "Points", "Point $i")
    set_attribute(gpe, "Indices", "$i" )
end


save_file(xmlSeries, joinpath(fishDir, "$(nNeurons)cell_sequential.xml"))
save_file(xmlGpl, joinpath(fishDir, "$(nNeurons)cell_sequential.gpl"))



## Then use image-block ripping, and we open 
# TODO replace with streaming...?


############# Analysis ##################
## load in results
# only triggered the first 120 trials, oops!!
# fishDir = "/media/tyler/tyler_benster/2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2"
# fishDir = "/data/dlab/b115/2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2"
tseriesName = "TSeries-1024cell-32concurrent-038"
# tseriesName = "TSeries-1024cell-32concurrent-043"
# tseriesName = "TSeries-cstoner_n64_b2_r8-043"
# tseriesName = "TSeries_1024cell_32concurrernt-045"

# 2021-01-11
tseriesName = "TSeries-1024cell-32concurrent-043"

tseriesDir = joinpath(fishDir, tseriesName)
slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

tseries = loadTseries(tseriesDir);
(H, W, Z, T) = size(tseries)
##
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]

# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
Z = size(etlVals,1)
volRate = frameRate / Z
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

# TODO: is stimStop wrong...? Seems too long....
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
allWithin(diff(stimStartIdx),0.05)
# adjust for expected number of stimuli....
@assert length(stimStartIdx) == 32
allWithin(stimEndIdx .- stimStartIdx,0.05)
avgStimDuration = mean(stimEndIdx .- stimStartIdx)
# adjust for volume rate / stim duration...
@assert 30 > avgStimDuration > 3


# Assumes no sequence of stim
# not be needed if in same session / avoid clobber
if ~(@isdefined target_groups)
    target_groups = [mat["cfg"]["maskS"]["targets"][1]
        for mat in matread.(findMatGroups(slmExpDir))]
end
    

nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
##
@assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}
is1024 = size(tseries,1)==1024
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
    is1024=is1024, zOffset=zOffset)

cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)

stimLocs = map(CartesianIndex ∘ Tuple, eachrow(cells[:,[2,1]]))
first(cells)

## emergency selection of period without motion...
@assert tif840path == "/mnt/deissero/users/tyler/b115/2020-12-17_h2b6s_chrmine-kv2.1_7dpf/fish2/SingleImage-840nm-1024-019/SingleImage-840nm-1024-019_Cycle00001_Ch3_000001.ome.tif"
lateImage = dropdims(mean(tseries[:,:,stimPlane,end-500:end], dims=3), dims=3)
start = 1000
earlyImage = dropdims(mean(tseries[:,:,stimPlane,start:start+500], dims=3), dims=3)
imshow(cat(reshape(lateImage,size(lateImage)..., 1),
     reshape(earlyImage,size(earlyImage)..., 1), dims=3))

##
stimPlane = Int(mean(map(x->mean(x[:,3]), targetsWithPlaneIndex)))
avgImage = dropdims(mean(tseries[:,:,stimPlane,:], dims=3), dims=3)
# avgImage = dropdims(mean(tseries[:,:,stimPlane,1:1000], dims=3), dims=3)
# avgImage = dropdims(mean(tseries[:,:,stimPlane,end-1000:end], dims=3), dims=3)
# avgImage = dropdims(mean(tseries[:,:,stimPlane,end-1300:end], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
avgImageAdj = RGB.(avgImageAdj)
channelview(avgImageAdj)[[1,3],:,:,:] .= 0
@warn "hardcode for 512x512"
microscope_units = (2* 0.6299544139175637μm, 2* 0.6299544139175637μm, 2.0μm)
# assume that "7um" spiral galvo signal is calibrated for 16x, not 25x
# "16x" is actually 14.4x
targetSizePx = 7μm * (14.4/25) / microscope_units[1]
avgImgWithTargets = addTargetsToImage(copy(avgImageAdj), cartIdx2Array(stimLocs),
    targetSize=targetSizePx)
save(joinpath(plotDir,"$(tseriesName)_avgImgWithTargets.png"), avgImgWithTargets)
##
cellDf_f = []
# @assert size(stimStartIdx,1)==size(stimLocs,1) # only did 120 by accident
winSize = Int(ceil(3*volRate))
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

sum(cellDf_f.>0.1), sum(cellDf_f.>0.4), sum(cellDf_f.>1.)
# imshow(addTargetsToImage(avgImage, stimLocs[[best],:]))
# imshow(tseries[:,:,1,stimStartIdx[idx]-winSize:stimEndIdx[idx]+winSize])
##
cmax = maximum(abs.(cellDf_f))
p_df_map = Gadfly.with_theme(:dark) do
    Gadfly.plot(cells[cells.df_f.<0.1,:], x=:x, y=:y, Gadfly.Geom.point, size=[targetSizePx], color=:df_f,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true, fixed=true))
end

img = PNG(joinpath(plotDir, "$(expName)_df_f_map_below10.png"), 6inch, 5inch)
Gadfly.draw(img, p_df_map)
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
function plotStim(tseries,roiMask,cells, idx, volRate; before=30, after=60)
    roiM = roiMask[idx]
    x, y, z, stimStart, stimStop = cells[idx,:]
    plotRange = stimStart-before:stimStop+after
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    fluorescentTrace = imageJkalmanFilter(fluorescentTrace)
    p = plot(plotRange/volRate, fluorescentTrace, left_margin=50px)
        # xlabel="time (s)", ylabel="fluorescence")
    plot!(Lensman.rectangle((stimStop-stimStart)/volRate,maximum(fluorescentTrace),
        stimStart/volRate,0), opacity=.5, label="stim")
    p
end
before = Int(ceil(volRate*10))
after = Int(ceil(volRate*20))
# best
plots = []
nplots = 64

# nplots = 8
for idx in rankings[end:-1:end-nplots+1]
    push!(plots, plotStim(tseries,roiMask,cells,idx, volRate, before=before, after=after))
end
p = plot(plots..., layout=(16, 4), legend=false, size=(1024,1024*2))
# p = plot(plots..., layout=(2, 4), legend=false, size=(1024,1024/4))

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
selectedStimLocs = stimLocs[rankings[end:-1:end-64+1]]
cells[rankings[end:-1:end-64+1],:]

## select cells with df/f >= 1.0
# min_df_F = 1.0
# maxIdx = searchsortedfirst(sort(cells, :df_f, rev=true).df_f, min_df_F, rev=true)
# selectedStimLocs = sample(stimLocs[rankings[end:-1:end-maxIdx+2]], 64, replace=false)

## horrible rescue for motion (please just don't)
validCellIdxs = cells.stimStart .>=1000
lateCells = cells[validCellIdxs,:]
cellsWeWantToStim = lateCells[sortperm(lateCells.df_f, rev=true),:][1:64,:]
selectedStimLocs = [CartesianIndex(row.y,row.x) for row in eachrow(cellsWeWantToStim)]

##
nCells = 64
base = 2
nReps = 2
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
count_concurrency(stonerGroups)
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

name = "cstoner_n64_b2_r2"
outname = joinpath(slmMaskOutputDir, name)
# save("$(outname).tif", im)

create_slm_stim(target_groups,
    joinpath(slmMaskOutputDir, name))

##
addTargetsToImage(copy(avgImageAdj), cartIdx2Array(selectedStimLocs),
    targetSize=targetSizePx)
