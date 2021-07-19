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

ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
else
    # tseriesRootDir = "/data/dlab/b115"
    tseriesRootDir = "/scratch/b115"
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

tseriesDir = "$tseriesRootDir/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"

avg_stim_h5_path = "$(tseriesDir)_avgStim_lstm.h5"
fish_dir = joinpath(splitpath(avg_stim_h5_path)[1:end-1]...)
plot_dir = joinpath(fish_dir, "plots-denoised")
if ~isdir(plot_dir)
    mkdir(plot_dir)
end
##
avgStim = h5read(avg_stim_h5_path, "/block1");
##
imshow(avgStim)
##
imshow(avgStim[:,:,6,16,:])
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
# @assert Z == 1

# cmax = 4
# cmin = -0.75
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
stimNum = 16

# f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
# TODO: should window be immediately after stim..? right now it's
# ~2 seconds after stim...?
f = mean(avgStim[:,:,:,:,end-window+1:end],dims=5)[:,:,:,:,1]
f0 = mean(avgStim[:,:,:,:,1:window],dims=5)[:,:,:,:,1]
df = f - f0
df_f = df./f0
nplots = 15
@assert nplots == nStimuli - 1

figW,figH = (figB*5, figB*3)
global fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(3,5, figsize=(figW,figH))
axs = permutedims(axs,(2,1))
z = 3
for s in 2:nStimuli
    ax = axs[s-1]
    global cim = ax.imshow(df_f[:,:,z,s], cmap="RdBu_r",
        norm=cnorm)
    n = s*2
    ax.set_title("$n cells")
    ax.set_axis_off()
end

# axs[1].imshow(stim_masks[:, :,z,stimNum], cmap="gray")

# cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
# plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
#     norm=cnorm)
# may need to adjust if colorbar is cutoff
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.91, 0.15, 0.0075, 0.7])
# cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
cbar = fig.colorbar(cim, cax=cbar_ax)
h = cbar.ax.set_ylabel("Δf/f",rotation=270)
# h.set_rotation = 180
fig

fig.savefig(joinpath(plot_dir,"titration_oneplane_df_f.svg"),
    dpi=600)
fig.savefig(joinpath(plot_dir,"titration_oneplane_df_f.png"),
    dpi=600)

##
imshow(df_f[:,:,16,3])
## Choose mask

img = Gray.(df_f[:,:,3,16])
img
##
draw!(img, Ellipse(
    Point(160,238), 35, 60))
img
##
function make_ellipse_mask(x,y,a,b,thesize=size(df_f)[1:2])
    mask = Gray.(zeros(Bool,thesize))
    draw!(mask, Ellipse(
        Point(x,y), a, b))
    reinterpret(Bool,mask)
end
medulla_mask = Gray.(zeros(Bool,size(df_f)[1:2]))
draw!(medulla_mask, Ellipse(
    Point(160,238), 35, 60))
medulla_mask
#
medulla_mask = reinterpret(Bool,medulla_mask);
##
medulla_sum = zeros(16)
z = 3
for s in 1:nStimuli
    medulla_sum[s] = mean(df_f[:,:,z,s][medulla_mask])
end
##

img = Gray.(df_f[:,:,3,16])
img

# draw!(img, Ellipse(
#     Point(388,268), 25, 25))
# img

draw!(img, Ellipse(
    Point(388,188), 25, 25))
img
##
rhab_mask = make_ellipse_mask(388,268,25,25)
lhab_mask = make_ellipse_mask(388,188,25,25)
Gray.(rhab_mask)
##
hab_sum = zeros(16)
z = 3
n_pix = sum(lhab_mask) + sum(rhab_mask)
for s in 1:size(avgStim)[end-1]
    lhab_sum = sum(df_f[:,:,z,s][lhab_mask])
    rhab_sum = sum(df_f[:,:,z,s][rhab_mask])
    hab_sum[s] = (lhab_sum + rhab_sum)/n_pix
end
plot(1:2:32,hab_sum, label="habenula",legend=:topleft)
plot!(1:2:32,medulla_sum, label="medulla")
xlabel!("# of stimulated habenular neurons")
ylabel!("Δf/f")
savefig(joinpath(plot_dir,"titration_region_trace.svg"))
savefig(joinpath(plot_dir,"titration_region_trace.png"))
##
masks = RGB.(zeros(Bool,size(df_f)[1:2]))
masks[medulla_mask] .= RGB(250/255,104/255,0)
masks[lhab_mask] .= RGB(109/255,135/255,100/255)
masks[rhab_mask] .= RGB(109/255,135/255,100/255)
masks

## look at early vs late

analysis_name = ""
# analysis_name = "kalman"

# tyh5Path = tseriesDir * ".ty.h5"
tyh5Path = tseriesDir * "2021-06-21_6pm.ty.h5"

if occursin("freq", tseriesDir)
    exp_param = :stimFreq
elseif occursin("power", tseriesDir)
    exp_param = :laserPower
end

if tseriesDir[end] == "/"
    tseriesDir = tseriesDir[1:end-1]
end


if occursin("freq", tseriesDir)
    exp_param = :stimFreq
elseif occursin("power", tseriesDir)
    exp_param = :laserPower
end

if tseriesDir[end] == "/"
    tseriesDir = tseriesDir[1:end-1]
end

fishDir = joinpath(splitpath(tseriesDir)[1:end-1]...)
expName = splitpath(tseriesDir)[end]
recording_folder = splitpath(tseriesDir)[end-2]
fish_name = splitpath(tseriesDir)[end-1]
tylerSLMDir = joinpath(fishDir, "slm")



@info "using ty.h5 file"
# tseries = read_tyh5(tyh5Path)
h5 = h5open(tyh5Path,"r")
dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm"
# dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden_init_from_pretrained-2021-06-21_6pm"
h5, tseries = lazy_read_tyh5(tyh5Path, dset);
plotDir = joinpath(fishDir, "plots-denoised")
avgStimStr = "_avgStim_lstm.h5"

##
voltageFiles = glob("*VoltageRecording*.csv", tseriesDir)

# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
Z = size(etlVals,1)
volRate = frameRate / Z
try
    global trialOrder, slmExpDir
    slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
    trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)
catch
    # not sure if needed
    # hack in case matlab code uses date program started not experiment run
    slmExpDir = joinpath(slmDir,Dates.format(expDate - Day(1), "dd-u-Y"))
    @info "try previous day"
    trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)
end

##
stimStartIdx, stimEndIdx, frameStartIdx = getStimTimesFromVoltages(voltageFiles[1], Z);
##
stim_of_interest = 16
mostCellStimIdxs = findall(trialOrder.==stim_of_interest)
nStimsPer = length(mostCellStimIdxs)
earlyStimStartIdx = stimStartIdx[mostCellStimIdxs[1:Int(nStimsPer/2)]]
earlyStimEndIdx = stimEndIdx[mostCellStimIdxs[1:Int(nStimsPer/2)]]
lateStimStartIdx = stimStartIdx[mostCellStimIdxs[Int(nStimsPer/2+1):end]]
lateStimEndIdx = stimEndIdx[mostCellStimIdxs[Int(nStimsPer/2+1):end]]


max_frames = stimStartIdx[1]-1 
max_time = (max_frames-1)/volRate # otherwise bad index issues

nseconds = 5
nseconds = minimum([nseconds, max_time])
@show "using $nseconds seconds."
pre = Int(ceil(nseconds*volRate))+1
post = Int(ceil(nseconds*volRate))+1

earlyStim = trialAverage(tseries, earlyStimStartIdx, earlyStimEndIdx,
    trialOrder[mostCellStimIdxs[1:Int(nStimsPer/2)]]; pre=pre, post=post);
##
lateStim = trialAverage(tseries, lateStimStartIdx, lateStimEndIdx,
    trialOrder[mostCellStimIdxs[1:Int(nStimsPer/2)]]; pre=pre, post=post);

##

function save_fig(avgStim, stimNum, prepend="")
    figB = 1.6

    window = minimum([Int(ceil(3*volRate)), max_frames])
    @assert (window <= post) & (window <= pre)
    cmax = 2.5
    cmin = -0.5
    # cmax = 0.5
    # cmin = -0.1
    # cmax = 4
    # cmin = -0.75
    # cmin = -4
    cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)

    f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    f0 = mean(avgStim[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    df = f - f0
    df_f = df./f0
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    if Z > 5
        figW,figH = (figB*Z/2, figB*2)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(2,Int(Z/2), figsize=(figW,figH))
    else
        figW,figH = (figB*Z, figB)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(1,Z, figsize=(figW,figH))
    end
    if Z==1
        axs = [axs]
    end
    # ax = axs[2]
    # ax = axs
    for z in 1:Z
        ax = axs[z]
        global cim = ax.imshow(df_f[:,:,z], cmap="RdBu_r",
            norm=cnorm)
        ax.set_axis_off()
        ax.set_title("$(Int(round(etlVals[z],digits=0)))μm")
        # this will make extra circles (1 extra per repetition...)
    end

    # may need to adjust if colorbar is cutoff
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.15, 0.0075, 0.7])
    # cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
    cbar = fig.colorbar(cim, cax=cbar_ax)
    path = joinpath(plotDir,"$(prepend)$(recording_folder)_$(fish_name)_$(expName)_$(analysis_name)_stim$stimNum")
    @show path*".svg"
    fig.savefig(path*".svg", dpi=600)
    fig.savefig(path*".png", dpi=600)
end
save_fig(earlyStim, 16, "early")
save_fig(lateStim, 16, "late")