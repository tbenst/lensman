ENV["DISPLAY"] = "localhost:11.0"
##
# using ImageView
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Unitful
import PyPlot
import Lensman: @pun, @assign
plt = PyPlot
matplotlib = plt.matplotlib
SMALL_SIZE = 7
MEDIUM_SIZE = 8
BIGGER_SIZE = 9
plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title
np = pyimport("numpy")
matscale = pyimport("matplotlib_scalebar.scalebar")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
init_workers()
##
resources = Resources();
r = Recordings[
    # "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
    "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
](;resources...
    # window_secs=10,
);
##
@pun (trial_average, vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, stim_end_idx, stim_start_idx, nstim_pulses, window_secs,
    fish_dir, tseries_units, lateral_unit, tseries_read_strategy
) = r;
@show nstim_pulses
@assert nstim_pulses == 10
tseries_read_strategy
##
imap = influence_map(trial_average, window_len);
##
targetSizePx = spiral_size(exp_date, lateral_unit)

if occursin("hwzt", string(tseries_read_strategy))
    analysis_name = "kalman_imap"
else
    error("meh")
end
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)
# window = window_len - 1
window = 6
# @assert Z == 1
cmax = 2.5
cmin = -0.5

cmax = 3
cmin = 0
# cmin = -0.2

ϵ = 0.0
# ϵ = 0.05

# cmax = 0.5
# cmin = -0.1
# cmax = 4
# cmin = -0.75
# cmin = -4
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)

ncol = 4
figB = Float64(uconvert(u"inch", 183u"mm") / (ncol*u"inch"))
figW,figH = (figB*ncol, figB*1.5)
fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(2, 4, figsize=(figW,figH))

ztotext(z) = "$(Int(round(etl_vals[z],digits=0)))"
zplanes = [5,10]
yrange = 40:420
felz_k = 50
felz_min = 10

for stimNum in 1:nStimuli
    df_f = imap[:,:,:,stimNum]
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    title_str = "{$(ztotext(zplanes[1])),$(ztotext(zplanes[2]))}μm"
    axs[1].set_title("targets z∈$(title_str)", fontsize=8)


    # show targets 
    if stimNum < 3
        ax = axs[stimNum]
    else
        ax = axs[1]
    end
    theSize = collect(size(df_f)[1:2])
    theSize[1] = length(yrange)
    ax.imshow(zeros(theSize...), cmap="RdBu_r", norm=cnorm)
    for (x,y,targetZ) in eachrow(unique(cells[cells.stimNum .== stimNum,[:x,:y,:z]]))
        if stimNum == 1
            z = 5
        elseif stimNum == 2
            z = 10
            text = "Raphe targets"
        end
        circle = matplotlib.patches.Circle((x,y-yrange[1]), targetSizePx, color="k",
            fill=false, lw=0.2)
        ax.add_patch(circle)
    end
    ax.set_axis_off()
    
    
    for (j,z) in enumerate(zplanes)
        ax_num = stimNum*2 + j
        ax = axs[ax_num]
        im = df_f[:,:,z]
        im = opening_median(im)
        segments = felzenszwalb(im, felz_k, felz_min)
        im[segments.image_indexmap .== 1] .= 0

        global cim = ax.imshow(im[yrange,:], cmap="viridis",
            clim=(cmin,cmax), interpolation="none")
        
        if ax_num == 3
            ax.set_title("2nd-order hab",fontsize=8)
        elseif ax_num == 5
            ax.set_title("raphe",fontsize=8)
        elseif ax_num == 7
            ax.set_title("control",fontsize=8)
        end


        ax.set_axis_off()
    end
end

um_per_px1 = tseries_units[2] / 1u"μm"
scalebar1 = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="upper right", font_properties=Dict("size" => 7))
axs[1].add_artist(scalebar1)

plt.tight_layout()
fig.subplots_adjust(right=0.9)
fig.subplots_adjust(wspace=0.01, hspace=0.01)
cbar_ax = fig.add_axes([0.91, 0.15, 0.0075, 0.7])
# cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
cbar = fig.colorbar(cim, cax=cbar_ax)
h = cbar.ax.set_ylabel("Δf/f",rotation=270, fontsize=7)
path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_8panel")
@show path*".svg"
fig.savefig(path*".svg", dpi=300)
fig.savefig(path*".png", dpi=300)
fig.savefig(path*".pdf", dpi=300)

fig
##
imshow(trial_average[:,:,5,:,:])

##
avg_tseries_path = joinpath(fish_dir, exp_name*"_kalman_trial_average.h5")
h5write(avg_tseries_path, "trial_average", trial_average);