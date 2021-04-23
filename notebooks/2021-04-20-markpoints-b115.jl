## sweep multiple stimulation freqencies, ie 5Hz - 30Hz
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
fishDir = "$tseriesroot/2021-04-20_h33r-chrmine_6dpf_6f/fish2"
# fishDir = "/scratch/b115/2021-03-09_h2b6s/fish2"
useRed = false

# greenpath = glob("*840*", fishDir)[end]
greenpath = glob("*920*", fishDir)[end]
# greenpath = glob("*920*", "/mnt/deissero/users/tyler/b115/2021-03-16_wt-chrmine_h2b6s/fish2")[end] # hack
greenpath = glob("*Ch3*.tif", greenpath)[end]
if useRed
    # redpath = glob("*1115*", fishDir)[1]
    redpath = glob("*1115*", fishDir)[end]
    redpath = glob("*Ch2*.tif", redpath)[end]
end

offset = float(uconvert(m, 0μm)) / m # no z offset for SLM2 as of 2021-02-15
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

neuron_locs = read_gpl(gpl_path, width=W, height=H, zoom=1)

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
imshow(addTargetsToImage(copy(rgb840), cartIdx2Array(neuron_locs),
    targetSize=targetSizePx))
##
"pan galvo for best stim quality with SLM."
function ideal_galvo_for_target(x, y)
    retX = 0.
    if x < 256
        retX = x + 50.
    else
        retX = x - 50.
    end
    retY = 0.
    if y < 256
        retY = y + 50.
    else
        retY = y - 50.
    end
    [retX, retY]
end


## k neurons per stim
k = 1
# we reverse to get (x,y)
targets_center_list = map(target->ideal_galvo_for_target(reverse(Tuple(target))...), neuron_locs)
if H == 1024
    targets_center_list ./= 2
end

target_groups = [vcat(cartIdxFunc.(locs, fill(offset, length(locs)))...)
    for locs in Iterators.partition(neuron_locs,k)]
@show size.(target_groups)

N = length(neuron_locs)
## Save files for SLM stim
powers = [1]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
name = "$(N)cell-$(k)concurrent-zoffset_$(zOffset)_$(nPowers)power"
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
    # targets_center_list=targets_center_list)

powerPerCell = 18
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 395mW)")


## trialOrder
# shuffle each block so we never have a stimuli more than once per block
trialOrder = vcat([randperm(length(target_groups)) for _ in 1:20]...)
write_trial_order(trialOrder, outname)