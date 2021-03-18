## sweep multiple stimulation freqencies, ie 5Hz - 30Hz
# ENV["DISPLAY"] = "localhost:13.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter, LightXML, LibExpat
import Gadfly
using Unitful: μm, m, s, mW
##
# fishDir = "/mnt/deissero/users/tyler/b115/2021-03-09_wt-chrmine-gc6f/fish1"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-03-15_wt-chrmine_h2b6s/fish1"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-03-16_h2b6s/fish1"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-03-16_wt-chrmine_h2b6s/fish2"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-03-16_rschrmine_h2b6s/fish3"
fishDir = "/mnt/deissero/users/tyler/b115/2021-03-16_h33r-chrmine_h2b6s/fish4"
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
im840 = ImageMagick.load(greenpath)
rgb840 = RGB.(imadjustintensity(adjust_gamma(im840, 2)))
channelview(rgb840)[[1,3],:,:] .= 0

if useRed
    im1000 = ImageMagick.load(redpath)
    rgb1000 = RGB.(imadjustintensity(adjust_gamma(im1000, 3.)))
    channelview(rgb1000)[[2,3],:,:] .= 0;
end
(H, W) = size(im840)
## read gpl targets
gpl_path = glob("*.gpl", fishDir)[1]

neuron_locs = read_gpl(gpl_path, width=W, height=H, zoom=1)

##

if size(im840,1)==512
    microscope_units = (2* 0.6299544139175637μm, 2* 0.6299544139175637μm, 2.0μm)
    cartIdxFunc = cartIdx2SeanTarget512
elseif size(im840,1)==1024
    microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
    cartIdxFunc = cartIdx2SeanTarget
else
    @error "unknown microscope units"
end

targetSizePx = 7μm * (14.4/25) / microscope_units[1]
imshow(addTargetsToImage(copy(rgb840), cartIdx2Array(neuron_locs),
    targetSize=targetSizePx))

## k neurons per stim
k = 8
target_groups = [vcat(cartIdxFunc.(locs, fill(offset, length(locs)))...)
    for locs in Iterators.partition(neuron_locs,k)]
@show size.(target_groups)

N = length(neuron_locs)
## Save files for SLM stim
powers = [1,0.5]
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
    # localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

powerPerCell = 33
println("Powers for power per cell of $powerPerCell: $(powerPerCell ./ 1000 .* 450mW)")