# need X11 forwarding for ssh, `ssh -X osprey`, then `echo $DISPLAY`, and fill below
ENV["DISPLAY"] = "localhost:11"
# ENV["DISPLAY"] = ":0"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, HDF5,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman
import Plots
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
using Base.Iterators: peel
import Unitful: μm
##
# can set in Julia extension settings for vscode, search "num threads"
stim1Dir = "/scratch/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab_raphe_40trial-045/"
# stim2Dir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab-raphe-20trial-part2-047/"

# VoltageRecording stopped short :/
# stim1Dir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2          .1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-039"

# stim1Dir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-lrhab_raphe_stim-40trial-038"
# next one TO-RUN
# stim2Dir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries_lrhab_raphe_40trial_part_2-040"

# stim1Dir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_lrhab_raphe_40trial-044"
# stim2Dir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_lrhab_raphe_40trial_part_2-046/"

# TO-RUN
# stim1Dir = "/oak/stanford/groups/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-039"
# tifDir = stim2Dir
tifDir = stim1Dir

# TODO: copy slmExpDir (see below) to, say,  `/scratch/SLM_files/`, so no dependency on network
slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"


expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
##
tseries = loadTseries(tifDir);
(H, W, Z, T) = size(tseries)
##
voltageFile = glob("*VoltageRecording*.csv", tifDir)[1]
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
@assert length(stimStartIdx) == 120

## single stim example
imshow(imadjustintensity(tseries[:,:,:,stimStartIdx[1]-1:stimEndIdx[1]+1]))

## read slm stim files
dataFolders = splitpath(tifDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
volRate = frameRate / Z

# if imaging many planes, may need to read another xml file since didn't get all planes in first file
@assert length(etlVals) == Z
@info "assume imaging from top-down"

#
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, _ = getTrialOrder(slmExpDir, expDate)


nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
@assert nTrials == size(stimStartIdx,1) # check TTL start-times match 

target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

# Array of matrices indicating x, y z (or was it y, x, z..?)
groupLocs = mapTargetGroupsToPlane(target_groups, etlVals)

# stimDuration = Int(ceil(2*volRate))
pre = Int(ceil(5*volRate))
post = Int(ceil(5*volRate))

## avg stim effect

# warmup JIT
# @time trialAverage(tseries, stimStartIdx[1:3], stimEndIdx[1:3], trialOrder[1:3],
#     pre=pre, post=post);
avgStim = trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder;
    pre=pre, post=post);
##
try
    h5write(joinpath(fishDir,expName*"_avgStim.h5"), "/block1", avgStim)
catch
end
##
avgStim = h5read(joinpath(fishDir,expName*"_avgStim.h5"), "/block1");
(H, W, Z, nStim) = size(avgStim)
## STA
window = Int(ceil(3*volRate))
@assert (window < post) & (window < pre)
cmax = 4
cmin = -0.75
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)

for stimNum in 1:3
    f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    f0 = mean(avgStim[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    df = f - f0
    global df_f = df./f0

    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    global fig = plt.figure(figsize=(50,10))
    plt.axis("off")
    # cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
    plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
        norm=cnorm)
    if ~isdir(joinpath(fishDir, "plots"))
        mkdir(joinpath(fishDir, "plots"))
    end
    fig.savefig(joinpath(fishDir, "plots",expName*"_stim$stimNum.png"))
end
fig
## use this plot for cropping colorbar
fig = plt.figure(figsize=(3,6), dpi=300)
plt.imshow(vcat([df_f[:,:,z] for z in 1:3]...), cmap="RdBu_r",
    norm=cnorm)
cbar = plt.colorbar()
cbar.ax.set_ylabel("df/f")
fig.savefig(joinpath(fishDir, "plots",expName*"_colorbar.png"))
fig

##
imshow(avgStim[:,:,:,1,:]) # trial type 1 (left hab)
imshow(avgStim[:,:,:,2,:]) # right hab
imshow(avgStim[:,:,:,3,:]) # raphe

## construct region masks
@info "hardcoded microscope units / targetSize"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
targetSizePx = (7μm * 14.4/25) / microscope_units[1]

# imshow(regionMasks)
regionMasks = constructGroupMasks(groupLocs, H, W, Z, targetSizePx=targetSizePx);

## clear stim tseries...
## does not free memory
tseries = zeros(1024,1024,5,2000)
tseries = nothing
GC.gc()
sleep(5)
GC.gc()

## Two cells: does free memory
tseries = zeros(1024,1024,5,2000)
tseries = nothing
GC.gc()
##
GC.gc()
##

## need to run cell twice for some reason (sleep not sufficient)
GC.gc()
##
# restingPreDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-pre-043"
# restingMidDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-mid-046/"
# restingPostDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-post-048"

# restingPreDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_resting_pre-042"
# restingMidDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_resting_mid-045"
# restingPostDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_resting_post-047"

restingPreDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSe"
restingMidDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries"
restingPostDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries"

restingPreDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-resting-pre-037"
# restingMidDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-resting-mid-039"
# restingPostDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries_resting_post-041"

# tifDir = restingMidDir
tifDir = restingPreDir
expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
tseries = loadTseries(tifDir)
(H, W, Z, T) = size(tseries)


## noise correlation
cmax = 0.1
lag = 0
for stimNum in 1:3
    global imageCorr = imageCorrWithMask(tseries[:,:,:,:],
        regionMasks[:,:,:,stimNum], lag=lag);
    # h5write(joinpath(fishDir,expName, "_imageCorr.h5"), "/imageCorr", imageCorr)
    global fig = plt.figure(figsize=(50,10))
    plt.axis("off")
    plt.imshow(hcat([imageCorr[:,:,z] for z in 1:Z]...), cmap="RdBu_r")
    plt.imshow(hcat([imageCorr[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
        clim=(-cmax, cmax))
    fig.savefig(joinpath(fishDir, "plots",expName*"_midCorr_$(lag)idx-lag_region$(stimNum)_cmax$(cmax).png"))
end
fig
## colorbar
fig = plt.figure(figsize=(4,6), dpi=300)
plt.imshow(vcat([tseries[:,:,z,1] for z in 1:3]...), cmap="RdBu_r",
    clim=(-cmax, cmax))
cbar = plt.colorbar()
cbar.ax.set_ylabel("correlation")
fig.savefig(joinpath(fishDir, "plots",expName*"_midCorr_colorbar.png"),
    bbox_inches="tight")
fig

##
trace = extractTrace(tseries, regionMasks[:,:,:,2]);

Plots.plot(trace)
##
# imshow(imageCorr)
