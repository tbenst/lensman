##
ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots
using Unitful: μm, m, s

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
##
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
    "/mnt/deissero/users/tyler/slm/masks/2020-11-30/fish1-256-cell")


## Then use image-block ripping, and we open 
# TODO replace with streaming...?



## load in results
# only triggered the first 120 trials, oops!!
tseriesDir = "/mnt/b115_mSLM/tyler/2020-11-30_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1/fov1-tectum/TSeries-64groups-4cell-stim-069"

H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tseriesDir)
addprocs(36)
@everywhere begin
    import Pkg
    Pkg.activate(".")
    using Lensman
end
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> Lensman._readTseriesTiffDir(S, framePlane2tiffPath, Z, T));

imshow(imadjustintensity(tseries[:,:,:,1:100]))
##
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
stimStartFrameIdx, stimEndFrameIdx = getStimTimesFromVoltages(voltageFile, Z)

## TODO:
# - calc df/f for all stimmed neurons
# - choose top 64
# - write out combo stim protocol
# targetH5path = glob("*_targets.h5", splitdir(tseriesDir)[1])[1]


targetH5path = "$(outname)_targets.h5"
targetsH5 = h5open(targetH5path, "r")
@read targetsH5 stimLocs

##
avgImage = dropdims(mean(tseries[:,:,1,1:1000], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
# addTargetsToImage(avgImageAdj, stimLocs)

##
# targetSize = 7μm/microscope_units[1]
@info "hardcoded microscope units"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
targetSize = 7μm/microscope_units[1]
roiDf_f = []
# @assert size(stimStartIdx,1)==size(stimLocs,1) # only did 120 by accident
winSize = 30
# TODO: don't hardcode frameRate (15)
delay=Int(ceil(15*4))
roiMask = []
for i = 1:size(stimStartFrameIdx,1)
    y,x = stimLocs[i,:]
    mask = Gray.(zeros(Bool, H,W))
    draw!(mask, Ellipse(CirclePointRadius(x,y,targetSize)))
    mask = findall(channelview(mask))
    push!(roiMask, mask)
    s, e = stimStartFrameIdx[i], stimEndFrameIdx[i]
    # hardcode 1 plane 
    f = mean(reinterpret(UInt16,tseries[mask, 1, e+delay:e+winSize]))
    f0 = mean(reinterpret(UInt16,tseries[mask, 1, s-winSize:s-1]))
    # df_f = (f-f0)/(f0+1e-8)
    df_f = (f-f0)
    push!(roiDf_f,df_f)
end

rankings = sortperm(roiDf_f)
best = rankings[end]
worst = rankings[1]
##
idx = rankings[end-10]
imshow(addTargetsToImage(avgImage, stimLocs[[best],:]))
imshow(tseries[:,:,1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize])

##
l = @layout [a; b]
idx = best
p1 = plot(mean(tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize],dims=1)', title="best")
for i in 1:5
    idx = rankings[end-i]
    plot!(mean(reinterpret(UInt16,tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize]),dims=1)')
end

idx = worst
p2 = plot(mean(reinterpret(UInt16,tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize]),dims=1)', title="worst")
for i in 1:5
    idx = rankings[1+i]
    plot!(mean(reinterpret(UInt16,tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize]),dims=1)')
end
plot(p1, p2, layout=l, legend=false)

# TODO:
# - need to add TrialOrder (that ought to really eff things up!!!)
# can pull out / refactor from 2020-11-30_analyze_3region
# - stimEndFrameIdx is VERY wrong... ( TTL based so it doesn't know when it ends...)
println( mean(stimEndFrameIdx .- stimStartFrameIdx))
##

neuron_loc