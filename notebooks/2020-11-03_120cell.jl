# unfinished, but has improved segmentation I think..?
##
true
##
import Pkg
Pkg.activate(".")
@warn raw"ensure using 25x AssignRoot in Z:\mSLM_B115\Main\multiSLM_Server\Startup, else SLM will be off"
@warn raw"ensure using Olympus 25x in Prairie View, else ETL will be off"
ENV["DISPLAY"] = "localhost:12"
##
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates, DataFrames, CSV, Statistics
    import Base.Threads.@threads
using Unitful: μm, m, s
sio = pyimport("scipy.io")
morph = pyimport("skimage.morphology")

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
tseriesFiles = joinpath.(tseriesDir, filter(x->(x[end-6:end]=="ome.tif") & occursin("Ch3", x),
    readdir(tseriesDir)))
tif0 = ImageMagick.load(tseriesFiles[1])
H, W = size(tif0)

tseries = zeros(Normed{UInt16,16}, H, W, size(tseriesFiles, 1))
@threads for t in 1:size(tseriesFiles,1)
    tseries[:,:,t] = ImageMagick.load(tseriesFiles[t])
end
# imshow(tseries[:,:,1:100])
##
slmDir = "/mnt/deissero/users/tyler/b115/2020-11-04_elavl3-chrmine-Kv2.1_h2b6s_7dpf/2020_11_5_0_37_51_"
trialOrder = CSV.File(open(read, joinpath(slmDir, "trialOrder.txt"))) |> DataFrame
trialOrder = convert(Array,trialOrder)[1,:]

##
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
voltages = CSV.File(open(read, voltageFile)) |> DataFrame
rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

##
plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
frameStarted = diff(voltages["frame starts"] .> std(voltages["frame starts"]).*3)
frameStartIdx = findall(frameStarted .== 1) .+ 1

##
@warn "hardcodes 5vol stim (1s @ 5Hz)"
stimStartIdx = findTTLStarts(voltages["respir"])[5:5:end]

## TODO:
# - calc df/f for all stimmed neurons
# - choose top 64
# - write out combo stim protocol
cartIdx2Array(neuron_locs)