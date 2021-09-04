ENV["DISPLAY"] = "localhost:11"
##
using ImageView
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed
import PyPlot
import Lensman: @pun, @assign
plt = PyPlot
matplotlib = plt.matplotlib
np = pyimport("numpy")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman

init_workers()
##
resources = Resources();
r = Recordings[
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...,
    # window_secs=10,
    tseries_read_strategy = :lazy_tiff,
    # tseries_read_strategy = :hwzt,
    # tseries_read_strategy = :lazy_hwzt,
    # tyh5_path="/data/dlab/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5",
    # tseries_dset=nothing
);
##
@pun (trial_average, vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, stim_end_idx, stim_start_idx, nstim_pulses, window_secs,
    tseries, fish_dir, tseries_read_strategy
) = r;
@show nstim_pulses
@assert nstim_pulses == 10
size(trial_average)
# rm(avgstim_path)
if tseries_read_strategy==:lazy_tyh5
    avgstim_path = joinpath(fish_dir,exp_name*"_lstm_avgStim.h5")
    h5write(avgstim_path, "/lstm", trial_average);
    analysis_name = "lstm_imap"
elseif tseries_read_strategy==:lazy_tiff
    avgstim_path = joinpath(fish_dir,exp_name*"_raw_avgStim.h5")
    h5write(avgstim_path, "/lstm", trial_average);
    analysis_name = "raw_imap"
elseif tseries_read_strategy==:lazy_hwzt
    avgstim_path = joinpath(fish_dir,exp_name*"_kalman_avgStim.h5")
    h5write(avgstim_path, "/kalman", trial_average);
    analysis_name = "kalman_imap"
end

trial_average = h5read(avgstim_path, "/kalman");
##
imap = influence_map(trial_average, window_len);
sum(imap, dims=[1,2,3])[1,1,1,:]
##
imshow(imap)
##
s = stim_start_idx[1]
imshow(shift_every_other_row(tseries[:,:,:,s-10:s+15],-3))
##
##
lateral_unit = microscope_lateral_unit(tseriesZ)
targetSizePx = spiral_size(exp_date, lateral_unit)
# analysis_name = "lstm"
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)
# window = window_len - 1
window = 6
# @assert Z == 1
# cmax = 2.5
# cmax = 5
# cmin = -0.5

cmax = 1.0
cmin = -0.2

ϵ = 0.0
# ϵ = 0.05

# cmax = 0.5
# cmin = -0.1
# cmax = 4
# cmin = -0.75
# cmin = -4
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)

@warn "df_f denominator ϵ may have changed"
# for stimNum in 1:nStimuli
for stimNum in 16:nStimuli
    df_f = imap[:,:,:,stimNum]
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    if tseriesZ > 5
        figW,figH = (figB*tseriesZ/2, figB*2)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(2,Int(tseriesZ/2), figsize=(figW,figH))
    else
        figW,figH = (figBtseriesZZ, figB)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(1,Z, figsize=(figW,figH))
    end
    if tseriesZ==1
        axs = [axs]
    end
    # ax = axs[2]
    # ax = axs
    for z in 1:tseriesZ
        ax = axs[z]
        global cim = ax.imshow(df_f[:,:,z], cmap="RdBu_r",
            norm=cnorm)
        ax.set_axis_off()
        ax.set_title("$(Int(round(etl_vals[z],digits=0)))μm")
        # this will make extra circles (1 extra per repetition...)
        for (x,y,targetZ) in eachrow(unique(cells[cells.stimNum .== stimNum,[:x,:y,:z]]))
            if z == targetZ
                circle = matplotlib.patches.Circle((x,y), targetSizePx, color="k",
                    fill=false, lw=0.4, alpha=0.3)
                ax.add_patch(circle)
            end
        end
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
    path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_stim$stimNum")
    @show path*".svg"
    # fig.savefig(path*".svg", dpi=600)
    # fig.savefig(path*".png", dpi=600)
    fig
end
fig


## Titration one-plane
figB = 1.6

nStimuli = size(trial_average)[end-1]
# nStimuli = size(long_avg)[end-1]
volRate = 3
max_frames = size(trial_average,5)
Z = size(trial_average,3)
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)

# cmax = 1.0
cmax = 2
cmin = -1.0
# cmin = -0.5
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

# f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
# TODO: should window be immediately after stim..? right now it's
# ~2 seconds after stim...?
f = mean(trial_average[:,:,:,:,end-window+1:end],dims=5)[:,:,:,:,1]
f0 = mean(trial_average[:,:,:,:,1:window],dims=5)[:,:,:,:,1]
df = f - f0
df_f = df./f0
nplots = 15
@assert nplots == nStimuli - 1

# figW,figH = (figB*4, figB*4)
figW,figH = (figB*4, figB*3.1)
global fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(4,4, figsize=(figW,figH))
axs = permutedims(axs,(2,1))
z = 3
for s in 1:nStimuli
    ax = axs[s]
    global cim = ax.imshow(df_f[50:end-100,1:end-50,z,s], cmap="RdBu_r",
        norm=cnorm)
    # global cim = ax.imshow(df_f[50:end-100,1:end-50,z,s], cmap="viridis", clim=(0.1,cmax))
    n = s*2
    # ax.set_title("$n cells")
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
plotpath = joinpath(plot_dir,"$(analysis_name)_rdbu_titration_oneplane_df_f")
fig.savefig("$(plotpath).svg",
    dpi=600)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=600)
@show plotpath*".png"
# fig
fig.subplots_adjust(wspace=0.01, hspace=0.01)
# plt.tight_layout()
fig