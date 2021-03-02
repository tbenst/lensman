## sweep multiple stimulation freqencies, ie 5Hz - 30Hz
# ENV["DISPLAY"] = "localhost:13.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter
import Gadfly
using Unitful: μm, m, s, mW
##
# fishDir = "/mnt/b115_data/tyler/2021-02-23_h2b6s/fish1/" # not run
# fishDir = "/mnt/b115_data/tyler/2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2"
fishDir = "/mnt/b115_data/tyler/2021-02-23_h2b6s/fish1_pt2"


tif840path = glob("*840*", fishDir)[end]
tif840path = glob("*Ch3*.tif", tif840path)[end]

offset = float(uconvert(m, 0μm)) / m # no z offset for SLM2 as of 2021-02-15
slmNum = 2
zOffset = offset * 1e6
expName = splitpath(tif840path)[end-1]
im840 = ImageMagick.load(tif840path)
rgb840 = RGB.(imadjustintensity(adjust_gamma(im840, 1.5)))
channelview(rgb840)[[1,3],:,:] .= 0
(H, W) = size(im840)
rgb840
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
# green only targetting...
candidateTargetLocs = greenCandidateTargetLocs
"$(length(greenCandidateTargetLocs)) green candidates"
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
    candArray = cartIdx2Array(candidateTargetLocs)
    
    # left eye or right eye
    eyeIdxs = (((candArray[:,1] .<= 287) .&
        (candArray[:,2] .>= 660))) .| (
        #     (candArray[:,1] .<= 342) .& (candArray[:,2] .>= 710)
        # )  .| (
            (candArray[:,1] .>= 717) .& (candArray[:,2] .>= 723)
        # )  .| (
        #     (candArray[:,1] .>= 771) .& (candArray[:,2] .>= 476)
        )
    candidateTargetLocs = candidateTargetLocs[(~).(eyeIdxs)]
end


imshow(addTargetsToImage(copy(rgb840), cartIdx2Array(candidateTargetLocs),
    targetSize=targetSizePx))

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
k = 4
target_groups = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
    for locs in Iterators.partition(neuron_locs,k)]

## buildup stims
# target_groups_1 = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
#     for locs in Iterators.partition(neuron_locs[1:16],1)]

# target_groups_2 = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
#     for locs in Iterators.partition(neuron_locs[1:32],2)]

# target_groups_4 = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
#     for locs in Iterators.partition(neuron_locs[1:64],4)]

# target_groups_8 = [vcat(cartIdx2SeanTarget.(locs, fill(offset, length(locs)))...)
#     for locs in Iterators.partition(neuron_locs[1:128],8)]

# target_groups = vcat(target_groups_1, target_groups_2, target_groups_4,
#     target_groups_8)

##
@show size.(target_groups)

N = length(neuron_locs)
## Save files for SLM stim
# nPower = 5
# name = "$(N)cell-$(k)concurrent-zoffset_$(zOffset)_$(nPower)powers"
# frequencies = [5,10,20,30]
powers = [1,0.75,0.5]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
# frequencies = [10,10]
# nFreq = length(frequencies)
# name = "$(N)cell-$(k)concurrent-zoffset_$(zOffset)_$(nFreq)freq"
name = "$(k)cell-concurrent_$(nPowers)powers_zoffset_$(zOffset)"
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
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 530mW)")


create_slm_stim(target_groups, outname,
    localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

# 32 groups, 1 stim every 5s = 1:40 per cycle. 8 reps x 3 powers = 64 min