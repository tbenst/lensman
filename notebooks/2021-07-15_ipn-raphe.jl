## compare 
# ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
using ImageView
import Gadfly
using Unitful: μm, m, s, mW
import Base.Threads.@threads
import PyPlot
import Base.size
plt = PyPlot
matplotlib = plt.matplotlib
Geom = Gadfly.Geom
Scale = Gadfly.Scale

ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
else
    tseriesRootDir = "/data/dlab/b115"
    # tseriesRootDir = "/scratch/b115"
    # tseriesRootDir = "/mnt/deissero/users/tyler/b115"
end

# newer
if ON_SHERLOCK
    slmDir = "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files"
else
    slmDir = "/mnt/deissero/users/tyler/b115/SLM_files"
    # slmDir = "/mnt/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment"
    # slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
    # older
    # slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"
    # slmDir = "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/"
end

tseriesDir = "$tseriesRootDir/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"

# avg_stim_h5_path = "$(tseriesDir)_avgStim_lstm.h5"
avg_stim_h5_path = "$(tseriesDir)_avgStim.h5"
fish_dir = joinpath(splitpath(avg_stim_h5_path)[1:end-1]...)
plot_dir = joinpath(fish_dir, "plots-denoised")
if ~isdir(plot_dir)
    mkdir(plot_dir)
end
##
avgStim = h5read(avg_stim_h5_path, "/block1");
##
imshow(avgStim[:,:,6,1,:])
##
figB = 1.6

nStimuli = size(avgStim)[end-1]
volRate = 3
max_frames = size(avgStim,5)
Z = size(avgStim,3)
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)

cmax = 1.0
cmin = -1.0
nseconds = 5
pre = Int(ceil(nseconds*volRate))+1
post = Int(ceil(nseconds*volRate))+1

window = minimum([Int(ceil(3*volRate)), max_frames])
@assert (window <= post) & (window <= pre)

# f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
# TODO: should window be immediately after stim..? right now it's
# ~2 seconds after stim...?
f = mean(avgStim[:,:,:,:,end-window+1:end],dims=5)[:,:,:,:,1]
f0 = mean(avgStim[:,:,:,:,1:window],dims=5)[:,:,:,:,1]
df = f - f0
df_f = df./(f0 .+ quantile(f0[:],0.2));

##
imshow(df_f) # 1 is IPN; 2 is raphe
size(df_f)
##
ipn_p = Gadfly.plot(x=df_f[:,:,:,1][:],
            Geom.histogram(bincount=50),
            Scale.y_log10)
raphe_p = Gadfly.plot(x=df_f[:,:,:,2][:],
            Geom.histogram(bincount=50),
            Scale.y_log10)

vstack(raphe_p, ipn_p)

##
imshow(df_f[:,:,6,2] .> 1.0)
imshow(df_f[:,:,6,1] .> 1.0)
##
img = RGB.(zeros(512,512,10))
# channelview(img)[1,:,:,:] .= (df_f[:,:,:,2] .> 1.0) .* 255
# channelview(img)[2,:,:,:] .= (df_f[:,:,:,1] .> 1.0) .* 255
channelview(img)[1,:,:,:] .= imadjustintensity(df_f[:,:,:,2])
channelview(img)[2,:,:,:] .= imadjustintensity(df_f[:,:,:,1])
imshow(img)

##
function make_ellipse_mask(x,y,a,b,thesize=size(df_f)[1:2])
    mask = Gray.(zeros(Bool,thesize))
    draw!(mask, Ellipse(
        Point(x,y), a, b))
    reinterpret(Bool,mask)
end

medulla_mask = Gray.(zeros(Bool,size(df_f)[1:2]))
draw!(medulla_mask, Ellipse(
    Point(235,238), 60, 60))
medulla_mask


##

mean(img,dims=3)[:,:,1]

##

"find overlap of IPN/raphe proj field in medulla."
function get_percent_overlap(threshold)
    mm = convert(Array{Bool}, medulla_mask)
    df = df_f .* mm
    raphe_hits = df[:,:,:,2] .> threshold
    ipn_hits = df[:,:,:,1] .> threshold
    sum(ipn_hits .& raphe_hits) / sum(ipn_hits .| raphe_hits)
end

"calc false discovery rate based on threshold compared to control"
function get_fdr(threshold)
    mm = convert(Array{Bool}, medulla_mask)
    df = df_f .* mm
    # df = df_f
    raphe_hits = sum(df[:,:,:,2] .> threshold)
    ipn_hits = sum(df[:,:,:,1] .> threshold)
    control_hits = sum(df[:,:,:,3] .> threshold)
    (control_hits)/(control_hits + (raphe_hits + ipn_hits)/2)
end
thresholds = 0.05:0.05:3.0
fdrs = map(get_fdr,thresholds)
Gadfly.plot(x=collect(thresholds), y=fdrs, Geom.line)
reasonable_thresh = collect(thresholds)[searchsortedfirst(fdrs,0.05,lt=(>))]

get_percent_overlap(reasonable_thresh)
##

mm3d = repeat(convert(Array{Bool}, medulla_mask),1,1,size(df,3))
nostim_df = selectdim(df_f,4,3)
medulla_nostim_px = nostim_df[mm3d]
reasonable_thresh = quantile(medulla_nostim_px,0.95) 5% false discovery rate 


map(get_percent_overlap, [collect(0.2:0.2:2.0)..., reasonable_thresh])

# @show sum(raphe_hits), sum(ipn_hits), sum(ipn_hits .| raphe_hits), sum(ipn_hits .& raphe_hits)

# reasonable_thresh = quantile(df_f[:,:,:,3][:],0.95)
# @show reasonable_thresh
# get_percent_overlap(reasonable_thresh)

## PLAN:
# choose medulla roi
# calculate df/f of:
    - medulla
    - raphe
    - IPN
    - medulla excluding targets

# what fraction of pixels in medulla show >0.2 df/f?
# what is overlap of raphe & IPN stimulation recruitment for recruited pixels?

plots