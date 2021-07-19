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
# fishDir = "$tseriesroot/2021-05-18_rsChRmine_h2b6s_6dpf/fish5"
# fishDir = "$tseriesroot/2021-05-18_hsChRmine_6f/fish2"
# fishDir = "$tseriesroot/2021-05-18_hsChRmine_h2b6s_6dpf/fish4"
# fishDir = "/scratch/b115/2021-03-09_h2b6s/fish2"

# fishDir = "$tseriesroot/2021-06-01_wt-chrmine_h2b6s/f1"
# fishDir = "$tseriesroot/2021-06-01_rsChRmine_h2b6s/fish5"
# fishDir = "$tseriesroot/2021-06-01_rsChRmine_h2b6s/fish3"
# fishDir = "$tseriesroot/2021-06-01_wt-chrmine_h2b6s/fish4"
# fishDir = "$tseriesroot/2021-06-02_rsChRmine-h2b6s/fish2"
# fishDir = "$tseriesroot/2021-06-02_wtChRmine_h2b6s/fish3"
# fishDir = "$tseriesroot/2021-06-08_rsChRmine_h2b6s/fish1" #bad
# fishDir = "$tseriesroot/2021-06-08_rsChRmine_h2b6s/fish2"
# fishDir = "$tseriesroot/2021-06-15_rsChRmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2021-06-29_hsChRmine_6f_6dpf/fish2"
# fishDir = "$tseriesroot/2021-07-14_rsChRmine_h2b6s_5dpf/fish1"
fishDir = "$tseriesroot/2021-07-14_rsChRmine_h2b6s_5dpf/fish2"

useRed = false

# greenpath = glob("*840*", fishDir)[end]
if length(glob("*840*", fishDir))>0
    greenpath = glob("*840*", fishDir)[end]
    # greenpath = glob("*840*", fishDir)[2]
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
# offset = float(uconvert(m, -30μm)) / m
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
println("Make sure to run `sudo chmod -R g+rw $fishDir`")
## read gpl targets
@assert length(glob("*.gpl", fishDir)) == 1 (@show glob("*.gpl", fishDir))
# gpl_path = glob("*.gpl", fishDir)[2]
gpl_path = glob("*.gpl", fishDir)[1]
# gpl_path = glob("*.gpl", fishDir)[3]

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
# name = "$(num_groups)groups_$(N)cells_$name_str"
# name = "round3_$(num_groups)groups_$(N)cells_$(name_str)"
expSLMdir = joinpath(fishDir,"slm")
outname = joinpath(expSLMdir, name_str)
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
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 280mW)")


## trialOrder
# shuffle each block so we never have a stimuli more than once per block
# trialOrder = vcat([randperm(length(target_groups)) for _ in 1:20]...)
ntransitions = 13
nattempts = 1000
trialOrder, successful = balanced_transition_order(num_groups, ntransitions,nattempts)
count_per_group = [sum(trialOrder.==n) for n=1:num_groups]
@assert all(abs.(count_per_group .- mean(count_per_group)) .<= 1)
if ~isfile(outname*"_trialOrder.txt")
    write_trial_order(trialOrder, outname)
else
    error("file already exists! refusing to clobber")
end
count_per_group
length(trialOrder)


## titration experiment
# lr_pairs = collect(zip(group_locs["right-hab"], group_locs["left-hab"]))
lr_pairs = collect(zip(group_locs["rhab"], group_locs["lhab"]))
# lr_pairs = lr_pairs[randperm(length(lr_pairs))]

target_groups = Matrix{Float64}[]

max_pairs = 16
titrations = collect(1:1:max_pairs)
# titrations = collect(2:2:max_pairs)

for i in titrations
    # collect 2*i neurons into a vector
    locs = vcat([collect(lr) for lr in lr_pairs[1:i]]...)
    push!(target_groups, 
        vcat(cartIdxFunc.(locs, fill(offset, length(locs)))...)
    )
end

@show size.(target_groups)

num_groups = length(target_groups)
## Save files for SLM stim
powers = [1]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
# name_str = "lrhab_titration"
name_str = "$(length(titrations))titrations_$(max_pairs*2)cells"
expSLMdir = joinpath(fishDir,"slm")
outname = joinpath(expSLMdir, name_str)
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

ntrials = 12
titration_trialOrder = vcat([randperm(num_groups) for _ in 1:ntrials]...)
write_trial_order(titration_trialOrder, outname)
@show length(titration_trialOrder)


## 




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

## Combo experiment....

nCells = 64
base = 2
nReps = 8
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
sameGroups = stimGroups
permGroups = stimGroups
stonerGroups = stimGroups
permutation = 1:nCells
for i in 0:nReps-1
    if i % 2 == 0
        # draw new permutation
        permutation = randperm(nCells)
    end
    if i % 2 == 1
        # choose new permutation that balances the previous
        stoner_perm = stonerPerm(nCells)
    else
        stoner_perm = collect(1:nCells)
    end
    stimGroupsPerm = map(g->perm(permutation,g),stimGroups)
    stimGroupsPerm = map(g->perm(stoner_perm,g),stimGroupsPerm)
    stonerGroups = vcat(stonerGroups, stimGroupsPerm)
end

stonerGroupsPerCell = calcGroupsPerCell(stonerGroups, nCells, base)


stonerGroups

stonerGroupsOpt, stonerGroupsPerCellOpt = randomSwaps(stonerGroups, stonerGroupsPerCell, calc_concurrency_score,
nCells, base, 10000)

count_concurrency(stonerGroupsOpt)
## GENERATE COMBO Experiment
all_cells = vcat(group_locs["rhab"], group_locs["lhab"])

##
k = Int(nCells / base)
target_groups = []
@warn "1024"
for group in eachrow(stonerGroups)
    groupStimLocs = all_cells[group]
    push!(target_groups, vcat(cartIdx2SeanTarget.(groupStimLocs, fill(offset, k))...))
end

## Save files for SLM stim
name_str = "cstoner_n$(nCells)_b$(base)_r$(nReps)"
outname = joinpath(expSLMdir, name_str)


create_slm_stim(target_groups, outname,
    # 9 is path "/scratch/" includes trailing
    localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

##
imshow(addTargetsToImage(copy(rgb840),
    cartIdx2Array(all_cells),
    targetSize=targetSizePx))