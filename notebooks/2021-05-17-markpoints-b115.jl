# now supports grouped markpoints
# ENV["DISPLAY"] = "localhost:13.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter, LightXML, LibExpat
import Gadfly
using Unitful: μm, m, s, mW
##
tseriesroot = "/scratch/b115"
# tseriesroot = "/mnt/deissero/users/tyler/b115"
# fishDir = "$tseriesroot/2021-03-09_wt-chrmine-gc6f/fish1"
# fishDir = "$tseriesroot/2021-03-15_wt-chrmine_h2b6s/fish1"
# fishDir = "$tseriesroot/2021-03-16_h2b6s/fish1"
# fishDir = "$tseriesroot/2021-03-16_wt-chrmine_h2b6s/fish2"
# fishDir = "$tseriesroot/2021-03-16_rschrmine_h2b6s/fish3"
# fishDir = "$tseriesroot/2021-03-16_h33r-chrmine_h2b6s/fish4"
# fishDir = "$tseriesroot/2021-03-30_wt-chrmine_6dpf_h2b6s/fish1"
# fishDir = "$tseriesroot/2021-04-19_wt-chrmine_5dpf_6f/fish1"
# fishDir = "$tseriesroot/2021-04-20_wt-chrmine_6dpf_6f/fish1"
# fishDir = "$tseriesroot/2021-04-20_h33r-chrmine_6dpf_6f/fish2"
# fishDir = "$tseriesroot/2021-05-11_wt-chrmine_6f_6dpf/fish1"
# fishDir = "$tseriesroot/2021-05-17_wt-chrmine_6f_6dpf/fish1"
# fishDir = "$tseriesroot/2021-05-18_rsChRmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2021-05-18_rsChRmine_h2b6s_6dpf/fish3"
fishDir = "$tseriesroot/2021-05-18_rsChRmine_h2b6s_6dpf/fish5"
# fishDir = "$tseriesroot/2021-05-18_hsChRmine_6f/fish2"
# fishDir = "$tseriesroot/2021-05-18_hsChRmine_h2b6s_6dpf/fish4"
# fishDir = "/scratch/b115/2021-03-09_h2b6s/fish2"
useRed = false

# greenpath = glob("*840*", fishDir)[end]
if length(glob("*840*", fishDir))>0
    # greenpath = glob("*840*", fishDir)[end]
    greenpath = glob("*840*", fishDir)[2]
else
    greenpath = glob("*920*", fishDir)[end]
end
# greenpath = glob("*920*", "/mnt/deissero/users/tyler/b115/2021-03-16_wt-chrmine_h2b6s/fish2")[end] # hack
greenpath = glob("*Ch3*.tif", greenpath)[end]
if useRed
    # redpath = glob("*1115*", fishDir)[1]
    redpath = glob("*1115*", fishDir)[end]
    redpath = glob("*Ch2*.tif", redpath)[end]
end

offset = float(uconvert(m, 0μm)) / m
@warn "using offset of $offset"
slmNum = 2
zOffset = offset * 1e6
expName = splitpath(greenpath)[end-1]
imgreen = ImageMagick.load(greenpath)
rgb840 = RGB.(imadjustintensity(adjust_gamma(imgreen, 2)))
channelview(rgb840)[[1,3],:,:] .= 0

if useRed
    im1000 = ImageMagick.load(redpath)
    rgb1000 = RGB.(imadjustintensity(adjust_gamma(im1000, 3.)))
    channelview(rgb1000)[[2,3],:,:] .= 0;
end
(H, W) = size(imgreen)
## read gpl targets
@assert length(glob("*.gpl", fishDir)) == 1 (@show glob("*.gpl", fishDir))
# gpl_path = glob("*.gpl", fishDir)[2]
gpl_path = glob("*.gpl", fishDir)[1]

group_locs = read_markpoint_groups(gpl_path, width=W, height=H, zoom=1)
# neuron_locs = read_gpl(gpl_path, width=W, height=H, zoom=1)

##

lateral_unit = microscope_lateral_unit(W)
microscope_units = (lateral_unit, lateral_unit, 2.0μm)
if W==512
    cartIdxFunc = cartIdx2SeanTarget512
elseif W==1024
    cartIdxFunc = cartIdx2SeanTarget
else
    @error "unknown microscope units"
end

targetSizePx = 7μm * (14.4/25) / microscope_units[1]
imshow(addTargetsToImage(copy(rgb840),
    cartIdx2Array(first(values(group_locs))),
    targetSize=targetSizePx))
##

target_groups = [vcat(cartIdxFunc.(locs, fill(offset, length(locs)))...)
    for locs in values(group_locs)]
@show size.(target_groups)

N = sum(length.(values(group_locs)))
num_groups = length(target_groups)
## Save files for SLM stim
powers = [1]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
name_str = join(keys(group_locs), "_")
name = "$(num_groups)groups_$(N)cells_$name_str"
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

create_slm_stim(target_groups, outname,
    # 9 is path "/scratch/" includes trailing
    localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

powerPerCell = 32
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 298mW)")


## trialOrder
# shuffle each block so we never have a stimuli more than once per block
# trialOrder = vcat([randperm(length(target_groups)) for _ in 1:20]...)
ntransitions = 10
nattempts = 1000
trialOrder, successful = balanced_transition_order(num_groups, ntransitions,nattempts)
count_per_group = [sum(trialOrder.==n) for n=1:num_groups]
@assert all(abs.(count_per_group .- mean(count_per_group)) .<= 1)
write_trial_order(trialOrder, outname)
length(trialOrder)

## one cell at a time

single_targets = [vcat(cartIdxFunc.([loc], fill(offset, 1))...)
    for loc in vcat(values(group_locs)...)]
idxs = vcat(collect(1:15), collect(33:47), collect(65:69))
# single_targets = single_targets[1:2:69] # only 5 control stims
single_targets = single_targets[idxs]
@show size.(single_targets)

N = length(single_targets)
num_groups = length(single_targets)
## Save files for SLM stim
powers = [1]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
name = "$(num_groups)groups_$(N)cells"
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

create_slm_stim(single_targets, outname,
    # 9 is path "/scratch/" includes trailing
    localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

powerPerCell = 32
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 298mW)")
##
ntrials = 20
single_trialOrder = vcat([randperm(N) for _ in 1:ntrials]...)
write_trial_order(single_trialOrder, outname)