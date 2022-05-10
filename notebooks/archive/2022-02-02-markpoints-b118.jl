# now supports grouped markpoints
# ENV["DISPLAY"] = "localhost:13.0"
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
fishDir = "$tseriesroot/2022-05-04_rschrmine_h2b6s_6dpf/fish1"
fishDir = "$tseriesroot/2022-05-04_rschrmine_h2b6s_6dpf/fish2"

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
@assert length(glob("*.gpl", fishDir)) == 1 (@show glob("*.gpl", fishDir))
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
xoffset, yoffset, estW, estH = getImagingROI(xml)
@assert estH ≈ H
@assert estW ≈ W

##
targetSizePx = 7μm * (14.4 / 25) / microscope_units[1]
adj_locs = vcat(
    cartIdx2Array(first(values(group_locs))),
    cartIdx2Array(collect(values(group_locs))[2])
)
adj_locs[:, 1] .-= Int(round(yoffset))
adj_locs[:, 2] .-= Int(round(xoffset))

# TODO: visualization does not work for ROI imaging?? points at wrong
imshow(addTargetsToImage(copy(rgb840),
    adj_locs,
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
expSLMdir = joinpath(fishDir, "slm")
outname = joinpath(expSLMdir, name_str)
if ~isdir(expSLMdir)
    mkdir(expSLMdir)
end
if isfile(outname * ".txt")
    error("file already exists! refusing to clobber")
end

# create_slm_stim(target_groups, outname,
#     localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
#     powers=collect(1:nPower)/nPower)

create_slm_stim(target_groups, outname,
    # 9 is path "/scratch/" includes trailing
    # localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\"),
    localToRemote=matpath -> "M:" * replace(matpath[15:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

powerPerCell = 32
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 280mW)")


## trialOrder
# shuffle each block so we never have a stimuli more than once per block
# trialOrder = vcat([randperm(length(target_groups)) for _ in 1:20]...)
ntransitions = 13
nattempts = 1000
trialOrder, successful = balanced_transition_order(num_groups, ntransitions, nattempts)
count_per_group = [sum(trialOrder .== n) for n = 1:num_groups]
@assert all(abs.(count_per_group .- mean(count_per_group)) .<= 1)
if ~isfile(outname * "_trialOrder.txt")
    write_trial_order(trialOrder, outname)
else
    error("file already exists! refusing to clobber")
end
count_per_group
length(trialOrder)
