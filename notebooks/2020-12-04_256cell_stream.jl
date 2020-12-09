##
ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose
import Gadfly
using Unitful: μm, m, s
##
# HOSTNAME = "171.65.101.131" # B113
HOSTNAME = "171.65.101.179" # B115
PORT = 1050
pv = connect(HOSTNAME, PORT)
println("connected")

# First, we read sizes of images to follow
bSize = read(pv, 8)
# reshape fails unless an Int64
C, Z, H, W = Int64.(reinterpret(Int16,bSize))
# TODO: right now only z=1 is functional...???
# at z=4, I get >50% more samples than expected...?
# perhaps "bad" pixels arent' being invalidated...?
# -> but no consistent visual effect if so...
# Z = 3
@show C, Z, H, W

currentVolume = Observable(zeros(Int16,C,H,W,Z))

numFrames = Observable(0)

images = []

# callback for new volumes
on(currentVolume) do vol
    numFrames[] += 1
    println("got volume $(numFrames[])")
    push!(images, vol)
end

# on(currentVolume) do vol
#     println("got volume with mean $(mean(vol))!")
# end

# start acquisition by sending TCP
tmpVol = zeros(Int16,C,H,W,Z)
println("writing to host")
write(pv, "Tseries")
# first image is from acquisition prior to now...
_ = read(pv, sizeof(Int16)*C*H*W)
# Next, we start reading images as they come in...
# TODO we're dropping frames, why?? PV python server
# says "acquired frame 1013"
# and we received 1012, but supposed to get 1024...
@async while isopen(pv)
    for z in 1:Z
        bArray = read(pv, sizeof(Int16)*C*H*W)
        println("read z $z")
        # we have to reassign currentVolume to trigger observable listener...
        arr = permutedims(reshape(reinterpret(Int16, bArray),C,W,H),[1,3,2])
        tmpVol[:,:,:,z] .= arr
    end
    currentVolume[] = copy(tmpVol)
end

nothing
##
close(pv) # warning: only run once done!
function zeroAdjust(image::Array{T,4}) where T<:Real
    # PMT offset
    # @assert sum(image.==0) == 0
    im = image .- 8192
    im[im.<0] .= 0
    # every other line, reverse column for roundtrip correction
    # for B115, need to start with odd line else fliplr
    im[:,1:2:end,:,:] .= im[:,1:2:end,end:-1:1,:]
    im
end

# compile
zeroAdjust(images[1])

avgImage = mean(zeroAdjust.(images))[2,:,:,1]
imshow(adjust_gamma(avgImage,0.5))
##
candNeuronCenterMask = findNeurons(avgImage,1.5,2,10);
candidateTargetLocs = findall(candNeuronCenterMask)
candidateTargets = copy(mask2ind(candNeuronCenterMask))

candidateTargets = 
# Gray.(candNeuronCenterMask)

img = RGB.(imadjustintensity(avgImage))
stim_points = zeros(Bool,size(img))
stim_points[candidateTargetLocs] .= true
stim_points = dilate(dilate(stim_points))
channelview(img)[[1,3],:,:,:] .= 0
channelview(img)[1,:,:,:] .= 0.25 * float(stim_points)
img

## Sample 128 neurons
nNeurons = 256
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
im

##
using Lensman
using Unitful: μm, m, s, uconvert

## to burn at etl=0 if using calibration circa 2020-11-16, need +30 offset
offset = float(uconvert(m, 30μm)) / m * 2
targets1 = [ 606. 636. offset;
           ] ./ 2

create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-30/fish1-single-cell")

##
# one cell at a time
target_groups = [cartIdx2SeanTarget(neuron_loc, offset)
    for neuron_loc in neuron_locs]

# k neurons per stim
k = 4
target_groups = [vcat(cartIdx2SeanTarget.(locs, fill(offset, k))...)
    for locs in Iterators.partition(neuron_locs,k)]


## Save files for SLM stim
outdir = "/mnt/deissero/users/tyler/slm/masks/2020-11-30/"
name = "fish1-256-cell"
outname = joinpath(outdir, name)
save("$(outname).tif", im)

HDF5.h5open("$(outname)_targets.h5", "w") do file
    stimLocs = copy(cartIdx2Array(neuron_locs))
    neuronLocs = copy(cartIdx2Array(candidateTargetLocs))
    HDF5.@write file neuronLocs
    HDF5.@write file stimLocs
    # maybe not needed/more confusing than anything...?
    HDF5.@write file targets
end

create_slm_stim(target_groups,
    outname)
    # "/mnt/deissero/users/tyler/slm/masks/2020-11-30/fish1-256-cell")


## Then use image-block ripping, and we open 
# TODO replace with streaming...?


############# Analysis ##################
## load in results
# only triggered the first 120 trials, oops!!
tseriesDir = "/mnt/b115_mSLM/tyler/2020-11-30_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1/fov1-tectum/TSeries-64groups-4cell-stim-069"
slmDir = "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/"
init_workers()
tseries = loadTseries(tseriesDir);
(H, W, Z, T) = size(tseries);
##
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
allWithin(diff(stimStartIdx),0.05)
@assert length(stimStartIdx) == 64
allWithin(stimEndIdx .- stimStartIdx,0.05)
avgStimDuration = mean(stimEndIdx .- stimStartIdx)
@assert 100 > avgStimDuration > 10
# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate = getExpData(xmlPath)
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

@warn "may not be needed if in same session / will clobber..."
# Assumes no sequence of stim
target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)

@assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}
cells = DataFrame(x=Int64[], y=Int64[], z=Int64[],
    stimStart=Int64[], stimStop=Int64[])
for (i,g) in enumerate(trialOrder)
    group = target_groups[g]
    start = stimStartIdx[i]
    endT = stimEndIdx[i]
    for (x,y,z) in eachrow(group)
        if W==1024
            x = x*2
            y = y*2
        end
        # hardcode z=1
        push!(cells, (x, y, 1, start, endT))
    end
end
first(cells)
##
stimLocs = map(CartesianIndex ∘ Tuple, eachrow(cells[:,[2,1]]))
##
# imshow(imadjustintensity(tseries[:,:,:,1:1000]))

## TODO:
# - calc df/f for all stimmed neurons
# - choose top 64
# - write out combo stim protocol


##
avgImage = dropdims(mean(tseries[:,:,1,1:1000], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
avgImageAdj = RGB.(avgImageAdj)
channelview(avgImageAdj)[[1,3],:,:,:] .= 0

addTargetsToImage(copy(avgImageAdj), cartIdx2Array(stimLocs))
##
# targetSize = 7μm/microscope_units[1]
@info "hardcoded microscope units"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
targetSize = 7μm/microscope_units[1]
# targetSize = 50
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
    draw!(view(mask,:,:,z), Ellipse(CirclePointRadius(x,y,targetSize)))
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
# addTargetsToImage(avgImageAdj, cartIdx2Array(stimLocs[rankings[end-31:end]]))
# addTargetsToImage(avgImageAdj, cartIdx2Array(stimLocs[rankings[end-31-32:end-32]]))

im = RGB.(avgImageAdj)
channelview(im)[[1,3],:,:,:] .= 0
im = addTargetsToImage(im, cartIdx2Array(stimLocs[rankings[end-63:end]]),
    channel="blue")
im = addTargetsToImage(im, cartIdx2Array(stimLocs[rankings[1:64]]),
    channel="red",alpha=0.6)
im

##
cmax = maximum(abs.(cellDf_f))
Gadfly.with_theme(:dark) do
    Gadfly.plot(cells, x=:x, y=:y, Gadfly.Geom.point, color=:df_f,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true))
end
##
function plotStim(tseries,roiMask,cells, idx)
    roiM = roiMask[idx]
    x, y, z, stimStart, stimStop = cells[idx,:]
    plotRange = stimStart-145:stimStop+145
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
plot(plots..., layout=(16, 4), legend=false, size=(870,1024*2))
## 33-64
plots = []
nplots = 32
for idx in rankings[end-32:-1:end-32-nplots+1]
    push!(plots, plotStim(tseries,roiMask,cells,idx))
end
plot(plots..., layout=nplots, legend=false, size=(1024,1024))

## 64-96
plots = []
nplots = 32
for idx in rankings[end-32:-1:end-32-nplots+1]
    push!(plots, plotStim(tseries,roiMask,cells,idx))
end
plot(plots..., layout=nplots, legend=false, size=(1024,1024))


## worst
plots = []
nplots = 32
for idx in rankings[1:nplots]
    push!(plots, plotStim(tseries,roiMask,cells,idx))
end
plot(plots..., layout=nplots, legend=false, size=(1024,1024))

## middle
plots = []
halfway = Int(size(cells,1)/2)-6+1
for idx in rankings[halfway:halfway+11]
    push!(plots, plotStim(tseries,roiMask,cells,idx))
end
l = @layout [a b c d; e f g h; i j k l]
plot(plots..., layout=l, legend=false)

##
plot(sum(tseries,dims=[1,2,3])[1,1,1,:])
##

Gadfly.plot(cells, x=:df_f, Gadfly.Geom.histogram)


### GENERATE COMBO Experiment
comboStimCells = sort(cells, :df_f, rev=true)[1:64,:]

k = 4
target_groups = [vcat(cartIdx2SeanTarget.(locs, fill(offset, k))...)
    for locs in Iterators.partition(neuron_locs,k)]


## Save files for SLM stim
outdir = "/mnt/deissero/users/tyler/slm/masks/2020-11-30/"
name = "fish1-256-cell"
outname = joinpath(outdir, name)
save("$(outname).tif", im)

HDF5.h5open("$(outname)_targets.h5", "w") do file
    stimLocs = copy(cartIdx2Array(neuron_locs))
    neuronLocs = copy(cartIdx2Array(candidateTargetLocs))
    HDF5.@write file neuronLocs
    HDF5.@write file stimLocs
    # maybe not needed/more confusing than anything...?
    HDF5.@write file targets
end

create_slm_stim(target_groups,
    "/mnt/deissero/users/tyler/slm/masks/2020-11-30/fish1-256-cell")

