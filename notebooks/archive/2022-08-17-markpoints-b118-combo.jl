# now supports grouped markpoints
ENV["DISPLAY"] = "localhost:12.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter, LightXML, LibExpat
import Gadfly
using Unitful: μm, m, s, mW
##
tseriesroot = "/mnt/b118_data/tyler"

# fishDir = "$tseriesroot/2022-02-02_rschrmine_h2b6s_7dpf/fish1"
# fishDir = "$tseriesroot/2022-02-09_rschrmine_h2b6s_7dpf/fish1"
# fishDir = "$tseriesroot/2022-03-08_rschrmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2022-03-11_rschrmine_h2b6s_9dpf/fish1"
# fishDir = "$tseriesroot/2022-03-11_rschrmine_h2b6s_9dpf/fish2"
# fishDir = "$tseriesroot/2022-03-31_rschrmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2022-03-31_rschrmine_h2b6s_6dpf/fish2"
# fishDir = "$tseriesroot/2022-05-04_rschrmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2022-05-19_rschrmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2022-06-07_wtchrmine_h2b6s_6dpf/fish1"
# fishDir = "$tseriesroot/2022-06-07_wtchrmine_h2b6s_6dpf/fish2"
fishDir = "$tseriesroot/2022-08-17_chrmine_h2b6s_6dpf/fish1"

useRed = false

if length(glob("*840*", fishDir)) > 0
    greenpath = glob("*840*", fishDir)[end]
    # greenpath = glob("*840*", fishDir)[2]
elseif length(glob("*860*", fishDir)) > 0
    greenpath = glob("*860*", fishDir)[end]
    # greenpath = glob("*860*", fishDir)[2]
else
    dirs = glob("*920*", fishDir)
    if length(dirs) == 1
        greenpath = dirs[end]
    else
        greenpath = dirs[end-1]
    end
end
greentifpath = glob("*Ch3*.tif", greenpath)[end]
@show greentifpath
if useRed
    # redpath = glob("*1115*", fishDir)[1]
    redpath = glob("*1115*", fishDir)[end]
    redpath = glob("*Ch2*.tif", redpath)[end]
end

# offset = float(uconvert(m, -30μm)) / m
offset = float(uconvert(m, 20μm)) / m
# offset = float(uconvert(m, 0μm)) / m
# offset = float(uconvert(m, 27μm)) / m
@warn "using offset of $offset"
slmNum = 1
zOffset = offset * 1e6
expName = splitpath(greentifpath)[end-1]
imgreen = ImageMagick.load(greentifpath)
rgb840 = RGB.(imadjustintensity(adjust_gamma(imgreen, 2)))
channelview(rgb840)[[1, 3], :, :] .= 0

if useRed
    im1000 = ImageMagick.load(redpath)
    rgb1000 = RGB.(imadjustintensity(adjust_gamma(im1000, 3.0)))
    channelview(rgb1000)[[2, 3], :, :] .= 0
end
(H, W) = size(imgreen)
H, W = 1024, 1024
@warn "FIXME: H,W needs to be the `base` sampling, e.g. 1024x1024, even if imaging ROI e.g. 1024x723"
# println("Make sure to run `sudo chmod -R g+rw $fishDir`")
## read gpl targets
# @assert length(glob("*.gpl", fishDir)) == 1 (@show glob("*.gpl", fishDir))
# gpl_path = glob("*.gpl", fishDir)[2]
gpl_path = glob("*.gpl", fishDir)[1]
# gpl_path = glob("*.gpl", fishDir)[5]
# gpl_path = glob("*.gpl", fishDir)[3]
@show gpl_path

group_locs = read_markpoint_groups(gpl_path, width=W, height=H, zoom=1)
# neuron_locs = read_gpl(gpl_path, width=W, height=H, zoom=1)

##
xmlfile = glob_one_file("*.xml", greenpath)
microscope_units = read_microns_per_pixel(xmlfile)

if W == 512
    cartIdxFunc = cartIdx2SeanTarget512
elseif W == 1024
    cartIdxFunc = cartIdx2SeanTarget
else
    @error "unknown microscope units"
end

@warn "wrong target size / not yet calibrated..?"

xml = parseXML(xmlfile)
xoffset, yoffset, estW, estH = Lensman.getImagingROI(xml)
@assert estH ≈ H
@assert estW ≈ W

##
targetSizePx = 7μm * (14.4 / 25) / microscope_units[1]
adj_locs = cartIdx2Array(first(values(group_locs)))
adj_locs[:, 1] .-= Int(round(yoffset))
adj_locs[:, 2] .-= Int(round(xoffset))

# TODO: visualization does not work for ROI imaging?? points at wrong
imshow(addTargetsToImage(copy(rgb840),
    adj_locs,
    targetSize=targetSizePx))

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
    stimGroupsPerm = map(g -> perm(permutation, g), stimGroups)
    stimGroupsPerm = map(g -> perm(stoner_perm, g), stimGroupsPerm)
    stonerGroups = vcat(stonerGroups, stimGroupsPerm)
end

stonerGroupsPerCell = calcGroupsPerCell(stonerGroups, nCells, base)


stonerGroups

stonerGroupsOpt, stonerGroupsPerCellOpt = randomSwaps(stonerGroups, stonerGroupsPerCell, calc_concurrency_score,
    nCells, base, 10000)

count_concurrency(stonerGroupsOpt)
## GENERATE COMBO Experiment
all_cells = vcat(group_locs["habenula"], group_locs["misc"])

##
k = Int(nCells / base)
target_groups = []
@warn "1024"
for group in eachrow(stonerGroups)
    groupStimLocs = all_cells[group]
    push!(target_groups, vcat(cartIdx2SeanTarget.(groupStimLocs, fill(offset, k))...))
end

## Save files for SLM stim
powers = [1]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
name_str = "cstoner_n$(nCells)_b$(base)_r$(nReps)"
expSLMdir = joinpath(fishDir, "slm")
outname = joinpath(expSLMdir, name_str)
if ~isdir(expSLMdir)
    mkdir(expSLMdir)
end


create_slm_stim(target_groups, outname,
    # 9 is path "/scratch/" includes trailing
    localToRemote=matpath -> "M:" * replace(matpath[15:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)
