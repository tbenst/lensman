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
resource = Resources()
##
# uri = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
uri = "2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-039"

r = Recordings[uri](
    ;resource...,
    tseries_read_strategy=:lazy_hwzt,
    tyh5_path="/data/dlab/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061_kalman.h5"
);
##
@pun (window_secs, trial_average, vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, nstim_pulses
) = r;
@show nstim_pulses
@assert nstim_pulses == 10
##
lateral_unit = microscope_lateral_unit(tseriesZ)
targetSizePx = spiral_size(exp_date, lateral_unit)

analysis_name = "raw_8-10s"
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)
# window = window_len - 1
window = 6
# @assert Z == 1
cmax = 2.5
cmin = -0.5

cmax = 1.0
cmin = -0.2

# ϵ = 0.0
ϵ = 0.05

# cmax = 0.5
# cmin = -0.1
# cmax = 4
# cmin = -0.75
# cmin = -4
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)

felz_k = 50
felz_min = 10

cmap = "inferno"
cmax = 1

@warn "df_f denominator ϵ may have changed"
for stimNum in 1:nStimuli
    f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    # f0 = mean(trial_average[:,:,:,stimNum,1:window-1],dims=4)[:,:,:,1]
    @warn "last 8-10s post stim.."
    f0 = mean(trial_average[:,:,:,stimNum,window_len-window:window_len],dims=4)[:,:,:,1]
    df = f - f0
    df_f = df./(f0 .+ ϵ)
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
        im = df_f[:,:,z]
        im = opening_median(im)
        segments = felzenszwalb(im, felz_k, felz_min)
        im[segments.image_indexmap .== 1] .= 0

        global cim = ax.imshow(im, cmap = cmap,
            # norm=cnorm,
            clim=(0,cmax),
            interpolation="none",)
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
end
fig