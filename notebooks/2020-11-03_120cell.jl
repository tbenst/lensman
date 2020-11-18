# unfinished, but has improved segmentation I think..?
##
true
##
@warn raw"ensure using 25x AssignRoot in Z:\mSLM_B115\Main\multiSLM_Server\Startup, else SLM will be off"
@warn raw"ensure using Olympus 25x in Prairie View, else ETL will be off"
ENV["DISPLAY"] = "localhost:11"

using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates, DataFrames, CSV, Statistics, Distributed, SharedArrays
    import Base.Threads.@threads
using Unitful: μm, m, s
sio = pyimport("scipy.io")
morph = pyimport("skimage.morphology")
##
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end
tifdir = "/mnt/deissero/users/tyler/b115/2020-11-04_elavl3-chrmine-Kv2.1_h2b6s_7dpf/SingleImage-840nm-010"

fishDir, _ = splitdir(tifdir)
tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") & occursin("Ch3", x),
    readdir(tifdir)))
plane = ImageMagick.load(tiff_files[1])
size(plane)

# TODO: make Olympus 25x in bruker have proper units....
@info "hardcoded microscope units"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)

# plane = AxisArray(plane, (:y, :x, :z), microscope_units)
adj_plane = adjust_histogram(imadjustintensity(plane), GammaCorrection(0.5))
# adj_plane = AxisArray(adj_plane, (:y, :x, :z), microscope_units)
@assert microscope_units[1] > 0.1μm
microscope_units

## neuron locations by plane
candNeuronCenterMask = findNeurons(plane,1.3)
candidateTargetLocs = findall(candNeuronCenterMask)
candidateTargets = copy(mask2ind(candNeuronCenterMask))

## Sample 128 neurons
nNeurons = 128
# nNeurons = 1
neuron_locs = sample(candidateTargetLocs, nNeurons, replace=false)


# all at once
targets = vcat(cartIdx2SeanTarget.(neuron_locs, zeros(Float64,nNeurons))...)
# one cell at a time
target_groups = [cartIdx2SeanTarget(neuron_loc, 0.)
    for neuron_loc in neuron_locs]

## Visualize stim targets
im = RGB.(adj_plane)
stim_points = zeros(Bool,size(im))
stim_points[neuron_locs] .= true
stim_points = dilate(dilate(stim_points))
channelview(im)[[1,3],:,:,:] .= 0
channelview(im)[1,:,:,:] .= 0.5*float(stim_points)
# imshow(im)
im

## Save files for SLM stim
outname = joinpath(fishDir, "120cell_sequential_1trial")
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
    "$(outname)_120targets")



## load in results
# only triggered the first 120 trials, oops!!
tseriesDir = "/mnt/deissero/users/tyler/b115/2020-11-04_elavl3-chrmine-Kv2.1_h2b6s_7dpf/TSeries-120seq-stim-042"

H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tseriesDir)
addprocs(18)
@everywhere begin
    import Pkg
    Pkg.activate(".")
    using Lensman
end
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> Lensman._readTseriesTiffDir(S, framePlane2tiffPath, Z, T));

imshow(imadjustintensity(tseries[:,:,:,1:100]))
##
slmDir = "/mnt/deissero/users/tyler/b115/2020-11-04_elavl3-chrmine-Kv2.1_h2b6s_7dpf/2020_11_5_0_37_51_"
trialOrder = CSV.File(open(read, joinpath(slmDir, "trialOrder.txt"))) |> DataFrame
trialOrder = convert(Array,trialOrder)[1,:]

voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
voltages = CSV.File(open(read, voltageFile)) |> DataFrame
rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

# plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
frameStarted = diff(voltages["frame starts"] .> std(voltages["frame starts"]).*3)
frameStartIdx = findall(frameStarted .== 1) .+ 1

@warn "hardcodes 5vol stim (1s @ 5Hz)"
stimDur = 5
stimStartIdx = findTTLStarts(voltages["respir"])[stimDur:stimDur:end]
stimEndIdx = findTTLEnds(voltages["respir"])[stimDur:stimDur:end]
stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimStartIdx]
stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimEndIdx]

## TODO:
# - calc df/f for all stimmed neurons
# - choose top 64
# - write out combo stim protocol
targetH5path = glob("*_targets.h5", splitdir(tseriesDir)[1])[1]
targetsH5 = h5open(targetH5path, "r")
@read targetsH5 stimLocs

##
avgImage = dropdims(mean(tseries[:,:,1,1:1000], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
# addTargetsToImage(avgImageAdj, stimLocs)

##
# targetSize = 7μm/microscope_units[1]
targetSize = 5μm/microscope_units[1]
roiDf_f = []
# @assert size(stimStartIdx,1)==size(stimLocs,1) # only did 120 by accident
winSize = 15
roiMask = []
for i = 1:size(stimStartIdx,1)
    y,x = stimLocs[i,:]
    mask = Gray.(zeros(Bool, H,W))
    draw!(mask, Ellipse(CirclePointRadius(x,y,targetSize)))
    mask = findall(channelview(mask))
    push!(roiMask, mask)
    s, e = stimStartFrameIdx[i], stimEndFrameIdx[i]
    # hardcode 1 plane 
    f = mean(reinterpret(UInt16,tseries[mask, 1, e+1:e+winSize]))
    f0 = mean(reinterpret(UInt16,tseries[mask, 1, s-winSize:s-1]))
    df_f = (f-f0)/(f0+1e-8)
    # df_f = (f-f0)
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
p1 = plot(mean(tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize],dims=1)')
for i in 1:5
    idx = rankings[end-i]
    plot!(mean(reinterpret(UInt16,tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize]),dims=1)')
end

idx = worst
p2 = plot(mean(reinterpret(UInt16,tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize]),dims=1)')
for i in 1:5
    idx = rankings[1+i]
    plot!(mean(reinterpret(UInt16,tseries[roiMask[idx],1,stimStartFrameIdx[idx]-winSize:stimEndFrameIdx[idx]+winSize]),dims=1)')
end
plot(p1, p2, layout=l, legend=false)