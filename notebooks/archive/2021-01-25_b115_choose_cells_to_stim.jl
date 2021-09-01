## characterize rsChRmine
ENV["DISPLAY"] = "localhost:12.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter
import Gadfly
using Unitful: μm, m, s, mW
##
# tif840path = "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish1/SingleImage-820nm-1024-018/SingleImage-820nm-1024-018_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish2/SingleImage-840nm-1024-021/SingleImage-840nm-1024-021_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/b115_data/tyler/2021-01-25_gcamp6f_6dpf/fish1/SingleImage-840nm-1024-023/SingleImage-840nm-1024-023_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/deissero/users/tyler/b115/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/SingleImage-840nm-1024-maxgdd-020/SingleImage-840nm-1024-maxgdd-020_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/b115_data/tyler/2021-01-25_chrmine_h2b6s_6dpf/fish1/SingleImage-840nm-1024-024/SingleImage-840nm-1024-024_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish3/SingleImage-840nm-1024-027/SingleImage-840nm-1024-027_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/b115_data/tyler/2021-01-26_rsChRmine_6f_7dpf/fish1/SingleImage-840nm-1024-023/SingleImage-840nm-1024-023_Cycle00001_Ch3_000001.ome.tif"
# tif840path = "/mnt/b115_data/tyler/2021-01-26_rsChRmine_6f_7dpf/fish2/SingleImage-840nm-1024-026/SingleImage-840nm-1024-026_Cycle00001_Ch3_000001.ome.tif"
tif840path = "/mnt/deissero/users/tyler/b115/2021-02-01_h33r_6f_5dpf/fish1/SingleImage-840nm-1024-023/SingleImage-840nm-1024-023_Cycle00001_Ch3_000001.ome.tif"
tif1000path = "/mnt/deissero/users/tyler/b115/2021-02-01_h33r_6f_5dpf/fish1/SingleImage-1000nm-1024-024/SingleImage-1000nm-1024-024_Cycle00001_Ch2_000001.ome.tif"
# tif920path = "/mnt/b115_data/tyler/2021-01-26_rsChRmine_6f_7dpf/fish1/SingleImage-920nm-1024-025/SingleImage-920nm-1024-025_Cycle00001_Ch3_000001.ome.tif"
# tif1000path = "/mnt/b115_data/tyler/2021-01-26_rsChRmine_6f_7dpf/fish1/SingleImage-1000nm-1024-024/SingleImage-1000nm-1024-024_Cycle00001_Ch2_000001.ome.tif"
# to burn at etl=0 if using calibration circa 2020-12-15, need +45 offset
# offset = float(uconvert(m, 45μm)) / m # prior to 2021 / starting on ...12/15...? should check...
offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
zOffset = offset * 1e6
fishDir = joinpath(splitpath(tif840path)[1:end-2]...)
expName = splitpath(tif840path)[end-1]
im840 = ImageMagick.load(tif840path)
rgb840 = RGB.(imadjustintensity(adjust_gamma(im840, 1.5)))
channelview(rgb840)[[1,3],:,:] .= 0
(H, W) = size(im840)
##
# im920 = ImageMagick.load(tif920path)
im1000 = ImageMagick.load(tif1000path)
channelview(im840)[[1,3],:,:] .= 0
rgb1000 = RGB.(imadjustintensity(adjust_gamma(im1000, 3.)))
channelview(im1000)[[2,3],:,:] .= 0
# rgb920 = RGB.(imadjustintensity(adjust_gamma(im920, 2.)))
# channelview(im920)[[1,3],:,:] .= 0
rgbBoth = rgb840
channelview(rgbBoth)[1,:,:] .= channelview(rgb1000)[1,:,:]
rgbBoth
## may want to now skip to Analysis section if not analyzing live
# thresh of 1.5 is reasonable for gcamp
# candNeuronCenterMask = findNeurons(im1000,thresh_adjust=1.5, featSize=4,maxiSize=8);
# thresh of 2.5 is reasonable for gcamp
redCandNeuronCenterMask = findNeurons(im1000,thresh_adjust=2.5, featSize=4,maxiSize=8);
redCandidateTargetLocs = findall(redCandNeuronCenterMask)
redEdgeIdxs = minimum(cartIdx2Array(redCandidateTargetLocs),dims=2) .<= 5
redEdgeIdxs = redEdgeIdxs[:,1] .| (maximum(cartIdx2Array(redCandidateTargetLocs),dims=2) .>= H-5)[:,1]
println("remove $(sum(redEdgeIdxs)) red edge indices")
redCandidateTargetLocs = redCandidateTargetLocs[(~).(redEdgeIdxs)]
##
greenCandNeuronCenterMask = findNeurons(im840,thresh_adjust=1.5, featSize=4,maxiSize=8);
greenCandidateTargetLocs = findall(greenCandNeuronCenterMask)
greenEdgeIdxs = minimum(cartIdx2Array(greenCandidateTargetLocs),dims=2) .<= 5
greenEdgeIdxs = greenEdgeIdxs[:,1] .| (maximum(cartIdx2Array(greenCandidateTargetLocs),dims=2) .>= H-5)[:,1]
println("remove $(sum(greenEdgeIdxs)) green edge indices")
greenCandidateTargetLocs = greenCandidateTargetLocs[(~).(greenEdgeIdxs)]
# green only targetting...
# candidateTargetLocs = greenCandidateTargetLocs
##
green_points = zeros(Bool,size(im840))
green_points[greenCandidateTargetLocs] .= true
greenTargets = findall(dilate(dilate(dilate(dilate(green_points)))))

candidateTargetLocs = intersect(greenTargets, redCandidateTargetLocs)
##

# if true
if false
    @warn "manual removal of eye idxs"
    candArray = cartIdx2Array(candidateTargetLocs)
    eyeIdxs = ((candArray[:,1] .<= 253) .&
        (candArray[:,2] .>= 669))
    candidateTargetLocs = candidateTargetLocs[(~).(eyeIdxs)]
end


# candidateTargets = copy(mask2ind(candNeuronCenterMask))

# addTargetsToImage(im840, candidateTargetLocs)

##
# img = copy(rgbBoth)
img = copy(rgb840)
stim_points = zeros(Bool,size(img))
stim_points[candidateTargetLocs] .= true
stim_points = dilate(dilate(stim_points))
# channelview(img)[[1,3],:,:,:] .= 0
# channelview(img)[3,:,:,:] .= 0.5 * float(stim_points)
channelview(img)[1,:,:,:] .= 0.3 * float(stim_points)
imshow(img)
println("found $(size(candidateTargetLocs,1)) potential targets")
## Sample 1024 neurons
nNeurons = 1024
# nNeurons = 1
if length(candidateTargetLocs) > nNeurons
    neuron_locs = sample(candidateTargetLocs, nNeurons, replace=false)     
else
    neuron_locs = candidateTargetLocs
end

# Visualize stim targets
img = RGB.(imadjustintensity(im840))
stim_points = zeros(Bool,size(img))
stim_points[neuron_locs] .= true
stim_points = dilate(dilate(stim_points))
channelview(img)[[1,3],:,:,:] .= 0
channelview(img)[1,:,:,:] .= 0.5*float(stim_points)
imshow(img)
img;


# one cell at a time
# target_groups = [cartIdx2SeanTarget(neuron_loc, offset)
#     for neuron_loc in neuron_locs]

# k neurons per stim
k = 32
target_groups = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
    for locs in Iterators.partition(neuron_locs,k)]
@show size.(target_groups)

N = length(neuron_locs)
## Save files for SLM stim
nPower = 4
name = "$(N)cell-$(k)concurrent-zoffset_$(zOffset)_$(nPower)powers"
outname = joinpath(fishDir, name)
if isfile(outname*".txt")
    error("file already exists! refusing to clobber")
end

create_slm_stim(target_groups, outname,
    localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
    powers=collect(1:nPower)/nPower)

24 .* collect(1:nPower)/nPower ./ 1000 .* 840mW

###### ROUND 2 ##########
## load TSeries & analyze cells.
# tyh5Path = "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish2/TSeries-1024cell-32concurrent-4power-043.ty.h5"
# tyh5Path = "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish3/TSeries-1024cell-32concurrent-4power-046.ty.h5"
# tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-25_gcamp6f_6dpf/fish1/TSeries-1024cell-32concurrent-4power-044.ty.h5"
tyh5Path = "/data/dlab/b115/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-25_gcamp6f_6dpf/fish1/TSeries-1024cell-32concurrent-4power-044.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043.ty.h5"
fishDir = joinpath(splitpath(tyh5Path)[1:end-1]...)
expName = replace(splitpath(tyh5Path)[end], ".ty.h5" => "")
tseriesDir = joinpath(fishDir, expName)

tseries = h5read(tyh5Path, "/imaging/raw")
@assert size(tseries,4)==1
tseries = permutedims(tseries, (2,1,3,4,5))
tseries = tseries[:,:,:,1,:];

(H, W, Z, T) = size(tseries)
@show (H, W, Z, T)
# look for motion artifact
avgImagingFrame = mean(tseries[:,:,:,1:10:end], dims=4)[:,:,:,1]
imshow(avgImagingFrame)
##


slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

# TODO: wrap more of cell in function to reduce inline code..?
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]

# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
numETLvals = size(etlVals,1)
@assert numETLvals == Z
volRate = frameRate / Z
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

# read power
@assert length(glob("*.txt", fishDir)) == 1 # if not, need to be careful to choose
slmTxtFile = glob("*.txt", fishDir)[1]
stimGroupDF = CSV.File(open(read, slmTxtFile), header=["filepath", "powerFraction"]) |> DataFrame
stimGroupDF = stimGroupDF[trialOrder,:]

stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
allWithin(diff(stimStartIdx),0.05)
# adjust for expected number of stimuli....
# @assert length(stimStartIdx) == 32
@show length(stimStartIdx)
allWithin(stimEndIdx .- stimStartIdx,0.05)
avgStimDuration = mean(stimEndIdx .- stimStartIdx)
# adjust for volume rate / stim duration...
@assert 30 > avgStimDuration > 3

# Assumes no sequence of stim
firstTargetGroup = matread.(findMatGroups(slmExpDir)[1])
lastTargetGroup = matread.(findMatGroups(slmExpDir)[end])
target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
# @assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}
is1024 = size(tseries,1)==1024
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
    is1024=is1024, zOffset=zOffset)

cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
@warn "hardcoded laser power"
slm1MaxPower = firstTargetGroup["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"]
##
slm1Power = 850mW
slm1powerPerCell = slm1Power * slm1MaxPower / 1000
cells[!, :slmNum] .= 1

cells[!, :laserPower] = round.(typeof(1.0mW), map(g->stimGroupDF.powerFraction[g], cells.stimGroup) .* slm1powerPerCell, digits=1)

stimLocs = map(CartesianIndex ∘ Tuple, eachrow(cells[:,[2,1]]))

stimPlane = Int(mean(map(x->mean(x[:,3]), targetsWithPlaneIndex)))
im840 = dropdims(mean(tseries[:,:,stimPlane,:], dims=3), dims=3)
im840Adj = adjust_gamma(imadjustintensity(im840), 0.5)
im840Adj = RGB.(im840Adj)
channelview(im840Adj)[[1,3],:,:,:] .= 0

if size(tseries,1)==512
    microscope_units = (2* 0.6299544139175637μm, 2* 0.6299544139175637μm, 2.0μm)
elseif size(tseries,1)==1024
    microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
else
    @error "unknown microscope units"
end

# assume that "7um" spiral galvo signal is calibrated for 16x, not 25x
# "16x" is actually 14.4x
targetSizePx = 7μm * (14.4/25) / microscope_units[1]
avgImgWithTargets = addTargetsToImage(copy(im840Adj), cartIdx2Array(stimLocs),
    targetSize=targetSizePx)
save(joinpath(plotDir,"$(expName)_avgImgWithTargets.png"), avgImgWithTargets)

winSize = Int(ceil(2*volRate))
# delay=Int(ceil(volRate*2))
delay=0
##
roiMask = constructROImasks(cells, H, W, Z, targetSizePx=targetSizePx+2);

cellsDF = makeCellsDF(tseries, cells, roiMask, winSize=winSize, delay=delay);
# meanDf_f = combine(groupby(cellsDF[cellsDF.laserPower .> 12mW,:], :cellID), :df_f => mean).df_f_mean
meanDf_f = combine(groupby(cellsDF, :cellID), :df_f => mean).df_f_mean
# maxDf_f = combine(groupby(cellsDF, :cellID), :df_f => maximum).df_f_maximum
@show [sum(cellsDF.df_f[cellsDF.laserPower .== p] .> 0.4)
    for p in sort(unique(cellsDF.laserPower))]
@show sum(meanDf_f .> 0.25)
cellIDrankings = sortperm(meanDf_f, lt=(>))
numNaN = sum(isnan.(cellsDF.df_f))
@assert numNaN == 0

@show sum(meanDf_f.>0.1), sum(meanDf_f.>0.2), sum(meanDf_f.>0.4), sum(meanDf_f.>1.)

for power in sort(unique(cellsDF.laserPower))
    println("$power: $(mean(cellsDF.df_f[cellsDF.laserPower .== power]))")
end
# only 26 are >40% for two stim in control...
# @show first(cellsDF[rankings[numNaN+1:end],:],50)

## viz what we actually stimulated...
targetSizePx = 7μm * (14.4/25) / microscope_units[1]
z = Int(mean(cellsDF.z))
im = RGB.(imadjustintensity(adjust_gamma(avgImagingFrame[:,:,z], 0.5)))
channelview(im)[[1,3],:,:] .= 0
avgImgWithTargets = addTargetsToImage(im, cartIdx2Array(stimLocs),
    targetSize=targetSizePx)

## viz STA (slow... memory depleted / swapping!)
error("memory pressure")
pre = Int(ceil(5*volRate))
post = Int(ceil(5*volRate))
avgStim = trialAverage(tseries, stimStartIdx, stimEndIdx,
    trialOrder; pre=pre, post=post);

imshow(avgStim[:,:,5,:,:])
##
Gadfly.plot(cellsDF, x="df_f",  ygroup="laserPower", Gadfly.Geom.subplot_grid(Gadfly.Geom.histogram))

##
before = Int(ceil(volRate*10))
after = Int(ceil(volRate*10))
# best
plots = []
# nplots = 40

nplots = 40
# for idx in rankings[end:-1:end-nplots+1]
# for idx in 1:nplots
@showprogress for cellID in cellIDrankings[1:nplots]
# for idx in randperm(nplots)
    indices = findall(cellsDF.cellID .== cellID)
    push!(plots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after))
end
# p = plot(plots..., layout=(16, 4), legend=false, size=(1024,1024*2))
p = plot(plots..., layout=(8,5),
    # size=(1024,512), legend=false)
    size=(1024*4,1024*2))

savefig(p, joinpath(plotDir,"$(expName)_top40_traces.png"))
savefig(p, joinpath(plotDir,"$(expName)_top40_traces.svg"))
p
##

chosenCellIdxs = cellIDrankings[[1, 2, 3, 6, 9, 10, 15, 17, 19, 21]]

chosenPlots = []
for cellID in chosenCellIdxs
    indices = findall(cellsDF.cellID .== cellID)
    push!(chosenPlots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after))
end
p = plot(chosenPlots..., layout=(2,5), size=(1024,512), legend=false)

## 

chosenTargetGroups = []
for (i,cellID) in enumerate(chosenCellIdxs)
    idx = findall(cellsDF.cellID .== cellID)[1]
    x, y = cellsDF[idx,[:x, :y]]
    push!(chosenTargetGroups, [x y offset])
end

# img = RGB.(imadjustintensity(im840))
# stim_points = zeros(Bool,size(img))
# stim_points[neuron_locs] .= true
# stim_points = dilate(dilate(stim_points))
# channelview(img)[[1,3],:,:,:] .= 0

# chosenTargetGroups

# Save files for SLM stim
name = "10cell-1concurrent-zoffset_$(zOffset)"
outname = joinpath(fishDir, name)
if isfile(outname*".txt")
    # error("file already exists! refusing to clobber")
end

create_slm_stim(chosenTargetGroups, outname,
    localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
    powers=[1])
