## sweep multiple stimulation freqencies, ie 5Hz - 30Hz
# ENV["DISPLAY"] = "localhost:13.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter
import Gadfly
using Unitful: μm, m, s, mW
##
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-15_wt_chrmine_gc6f/fish1"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-16_h2b6s_h33r_f0_6dpf"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-16_6f_h33r_f0_6dpf/fish2"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-16_h2b6s_wt-chrmine/fish3"
# fishDir = "/mnt/b115_data/tyler/2021-02-23_h2b6s/fish1/SingleImage-840nm-1024-027"
fishDir = "/mnt/deissero/users/tyler/b115/2021-03-02_h33r-6f-6dpf/fish1"
useRed = true

tif840path = glob("*840*", fishDir)[end]
tif840path = glob("*Ch3*.tif", tif840path)[end]
if useRed
    tif1000path = glob("*1000*", fishDir)[end]
    tif1000path = glob("*Ch2*.tif", tif1000path)[end]
end

offset = float(uconvert(m, 0μm)) / m # no z offset for SLM2 as of 2021-02-15
slmNum = 2
zOffset = offset * 1e6
expName = splitpath(tif840path)[end-1]
im840 = ImageMagick.load(tif840path)
rgb840 = RGB.(imadjustintensity(adjust_gamma(im840, 2)))
channelview(rgb840)[[1,3],:,:] .= 0

if useRed
    im1000 = ImageMagick.load(tif1000path)
    rgb1000 = RGB.(imadjustintensity(adjust_gamma(im1000, 3.)))
    channelview(rgb1000)[[1,3],:,:] .= 0;
end
(H, W) = size(im840)
##



## may want to now skip to Analysis section if not analyzing live
# thresh of 1.5 is reasonable for gcamp
# candNeuronCenterMask = findNeurons(im1000,thresh_adjust=1.5, featSize=4,maxiSize=8);
# thresh of 2.5 is reasonable for gcamp
##
greenCandNeuronCenterMask = findNeurons(im840,thresh_adjust=1.5, featSize=4,maxiSize=8);
greenCandidateTargetLocs = findall(greenCandNeuronCenterMask)
greenEdgeIdxs = minimum(cartIdx2Array(greenCandidateTargetLocs),dims=2) .<= 5
greenEdgeIdxs = greenEdgeIdxs[:,1] .| (maximum(cartIdx2Array(greenCandidateTargetLocs),dims=2) .>= H-5)[:,1]
println("remove $(sum(greenEdgeIdxs)) green edge indices")
greenCandidateTargetLocs = greenCandidateTargetLocs[(~).(greenEdgeIdxs)]


if useRed
    redCandNeuronCenterMask = findNeurons(im1000,thresh_adjust=2.5, featSize=4,maxiSize=8);
    redCandidateTargetLocs = findall(redCandNeuronCenterMask)
    redEdgeIdxs = minimum(cartIdx2Array(redCandidateTargetLocs),dims=2) .<= 5
    redEdgeIdxs = redEdgeIdxs[:,1] .| (maximum(cartIdx2Array(redCandidateTargetLocs),dims=2) .>= H-5)[:,1]
    println("remove $(sum(redEdgeIdxs)) red edge indices")
    redCandidateTargetLocs = redCandidateTargetLocs[(~).(redEdgeIdxs)]
    "$(length(redCandidateTargetLocs)) red candidates"
    red_points = zeros(Bool,size(im840))
    red_points[redCandidateTargetLocs] .= true
    redTargetCircs = findall(dilate(dilate(dilate(dilate(red_points)))))
    candidateTargetLocs = intersect(greenCandidateTargetLocs, redTargetCircs)
    println("$(length(candidateTargetLocs)) green&red candidates")
    # and sample more from green until we have 1024 targets
    extraCandidateTargetLocs = setdiff(greenCandidateTargetLocs, redTargetCircs)
    println("$(length(extraCandidateTargetLocs)) green-only candidates")
else
    # green only targetting...
    candidateTargetLocs = greenCandidateTargetLocs
    "$(length(greenCandidateTargetLocs)) green candidates"
end

if useRed
    redCandNeuronCenterMask = findNeurons(im1000,thresh_adjust=2.5, featSize=4,maxiSize=8);
    redCandidateTargetLocs = findall(redCandNeuronCenterMask)
    redEdgeIdxs = minimum(cartIdx2Array(redCandidateTargetLocs),dims=2) .<= 5
    redEdgeIdxs = redEdgeIdxs[:,1] .| (maximum(cartIdx2Array(redCandidateTargetLocs),dims=2) .>= H-5)[:,1]
    println("remove $(sum(redEdgeIdxs)) red edge indices")
    redCandidateTargetLocs = redCandidateTargetLocs[(~).(redEdgeIdxs)]
    "$(length(redCandidateTargetLocs)) red candidates"
end
##

if size(im840,1)==512
    microscope_units = (2* 0.6299544139175637μm, 2* 0.6299544139175637μm, 2.0μm)
elseif size(im840,1)==1024
    microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
else
    @error "unknown microscope units"
end

targetSizePx = 7μm * (14.4/25) / microscope_units[1]
imshow(addTargetsToImage(copy(rgb840), cartIdx2Array(candidateTargetLocs),
    targetSize=targetSizePx))
##


if true
# if false
    @warn "manual removal of eye idxs"
    eye_indices(arr) = (
            (arr[:,1] .<= 338) .& (arr[:,2] .>= 465)
        ) .| (
            (arr[:,1] .>= 642) .& (arr[:,2] .>= 450)
        )

    eyeIdxs = eye_indices(cartIdx2Array(candidateTargetLocs))
    candidateTargetLocs = candidateTargetLocs[(~).(eyeIdxs)]

    eyeIdxs = eye_indices(cartIdx2Array(extraCandidateTargetLocs))
    extraCandidateTargetLocs = extraCandidateTargetLocs[(~).(eyeIdxs)]
end


# imshow(addTargetsToImage(copy(rgb840), cartIdx2Array(candidateTargetLocs),
#     targetSize=targetSizePx))

# backupneuron_locs = copy(neuron_locs)
## Sample 1024 neurons
# nNeurons = 1024
nNeurons = 128
# nNeurons = 1
if length(candidateTargetLocs) >= nNeurons
    neuron_locs = sample(candidateTargetLocs, nNeurons, replace=false)
    # neuron_locs = sample(neuron_locs, nNeurons, replace=false)
else
    nToAdd = nNeurons - length(candidateTargetLocs)
    extras = sample(extraCandidateTargetLocs, nToAdd, replace=false)
    neuron_locs = vcat(candidateTargetLocs, extras)
end
##
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

## k neurons per stim
k = 8
target_groups = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
    for locs in Iterators.partition(neuron_locs,k)]
@show size.(target_groups)

N = length(neuron_locs)
## Save files for SLM stim
# nPower = 5
# name = "$(N)cell-$(k)concurrent-zoffset_$(zOffset)_$(nPower)powers"
# frequencies = [5,10,20,30]
frequencies = [5]
nFreq = length(frequencies)
powers = repeat([1], nFreq)
name = "$(N)cell-$(k)concurrent-zoffset_$(zOffset)_$(nFreq)freq"
expSLMdir = joinpath(fishDir,"slm")
outname = joinpath(expSLMdir, name)
if ~isdir(expSLMdir)
    mkdir(expSLMdir)
end
if isfile(outname*".txt")
    error("file already exists! refusing to clobber")
end

# create_slm_stim(target_groups, outname,
#     localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
#     powers=collect(1:nPower)/nPower)

powerPerCell = 28
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 300mW)")


create_slm_stim(target_groups, outname,
    # localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

