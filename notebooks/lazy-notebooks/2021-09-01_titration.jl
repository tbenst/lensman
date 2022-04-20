ENV["DISPLAY"] = "localhost:11"
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

matplotlib.rcParams["font.sans-serif"] = ["Arial", "sans-serif"]

init_workers()
resources = Resources();
##
r = Recordings[
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...,
    # LSTM tyh5 looks best
    # tseries_read_strategy = :lazy_hwzt,
    # tyh5_path="/data/dlab/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5",
    # tseries_dset=nothing
);
##
@pun (vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, stim_end_idx, stim_start_idx, nstim_pulses, window_secs,
    tseries, fish_dir, tseries_read_strategy, tseries_units, lateral_unit,
    region_masks_h5, zbrain_mask_names, imaging2zseries_plane, tseriesH, tseriesW,
    target_groups, zseries, tseries, tseriesT
) = r;
@show nstim_pulses
@assert nstim_pulses == 10
# rm(avgstim_path)

if occursin("tyh5", string(tseries_read_strategy))
    avgstim_path = joinpath(fish_dir,exp_name*"_lstm_avgStim.h5")
elseif occursin("tiff", string(tseries_read_strategy))
    avgstim_path = joinpath(fish_dir,exp_name*"_raw_avgStim.h5")
elseif occursin("hwzt", string(tseries_read_strategy))
    avgstim_path = joinpath(fish_dir,exp_name*"_kalman_avgStim.h5")
end


try 
    if occursin("tyh5", string(tseries_read_strategy))
        global analysis_name = "lstm_imap"
        global trial_average = h5read(avgstim_path, "/lstm");
    elseif occursin("hwzt", string(tseries_read_strategy))
        global analysis_name = "kalman_imap"
        global trial_average = h5read(avgstim_path, "/kalman");
    else
        global analysis_name = "raw_imap"
        # global trial_average = h5read(avgstim_path, "/lstm");
        global trial_average = h5read(avgstim_path, "/raw");
    end
catch
    for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println(stdout)
    end
    println("couldn't read $avgstim_path")
    error("oops")
    global trial_average = r[:trial_average]
    if tseries_read_strategy==:lazy_tyh5
        h5write(avgstim_path, "/lstm", trial_average);
    elseif tseries_read_strategy==:lazy_tiff
        h5write(avgstim_path, "/raw", trial_average);
    elseif tseries_read_strategy==:lazy_hwzt
        h5write(avgstim_path, "/kalman", trial_average);
    end
end
size(trial_average)
##
imap = influence_map(trial_average, window_len);
sum(imap, dims=[1,2,3])[1,1,1,:]
##
# imshow(imap)
##
s = stim_start_idx[1]
# imshow(shift_every_other_row(tseries[:,:,:,s-10:s+15],-3))
##
##
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

## brain region outlines
REGION_LIST = [
    # "Cerebellum"
    # "Rhombencephalon - Gad1b Cluster 1"
    # "Rhombencephalon - Gad1b Cluster 16"
    # "Rhombencephalon - Gad1b Cluster 2"
    # "Glyt2 Cluster 12"
    # "Glyt2 Stripe 3"
    # "Habenula"
    # "Tectum Stratum Periventriculare"
    # "Raphe - Superior"
    # "Subpallial Gad1b cluster"
    "Rhombencephalon -"
    "Mesencephalon -"
    # "Telencephalon -"
    # "Diencephalon -" # we will manually draw
]
region_masks = [L.read_first_mask(region_masks_h5, zbrain_mask_names,
imaging2zseries_plane, region; outline=false) for region in REGION_LIST];
# dilate on outside for thicker line

##
region_outlines = [dilate(L.mask2outline(L.unsalt_n_pepper(regm[2]; felz_min=5000)),[1,2])
    for regm in region_masks];
    # for regm in region_masks[[1]]];
outlines = cat(region_outlines...,dims=4);
outlines = maximum(outlines,dims=4)[:,:,:,1];
small_outlines = imresize(outlines, (tseriesH, tseriesW, tseriesZ));
Gray.(outlines[:,:,3])

#################################
## Titration one-plane
#################################
figB = 1.6

nStimuli = size(trial_average)[end-1]
# nStimuli = size(long_avg)[end-1]
volRate = 3
max_frames = size(trial_average,5)
Z = size(trial_average,3)
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)

cmax = 1.0
# cmax = 2
# cmin = -1.0
cmin = 0
# cmin = -0.5
nseconds = 5
pre = Int(ceil(nseconds*volRate))+1
post = Int(ceil(nseconds*volRate))+1

window = minimum([Int(ceil(3*volRate)), max_frames])
@assert (window <= post) & (window <= pre)
# @assert Z == 1

# cmax = 4
# cmin = -0.75
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
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
# two-column width (mm)
figB = Float64(uconvert(u"inch", (183-3)u"mm") / 4u"inch")
figW,figH = (figB*4, figB*3)
global fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(4,4, figsize=(figW,figH))
axs = permutedims(axs,(2,1))
z = 3
felz_k = 50
felz_min = 10
# show_outline = false
no_df = false
show_outline = true
# no_df = true
# yrange = 65:tseriesH-100
yrange = 85:tseriesH-120
# xrange = 60:tseriesW-90
xrange = 40:tseriesW-70
for s in 1:nStimuli
    ax = axs[s]
    im = df_f[:,:,z,s]
    im = opening_median(im)
    segments = felzenszwalb(im, felz_k, felz_min)
    im[segments.image_indexmap .== 1] .= 0
    im = im[yrange,xrange]
    # global cim = ax.imshow(df_f[50:end-100,1:end-50,z,s], cmap="RdBu_r",
    # global cim = ax.imshow(im, cmap="RdBu_r",
    if ~no_df 
        global cim = ax.imshow(im, cmap="viridis",
            clim=(cmin,cmax), interpolation="none")
            # norm=cnorm, interpolation="none")
    end
    # ax.imshow(outlines[:,:,z], alpha=outlines[:,:,z])
    # global cim = ax.imshow(df_f[50:end-100,1:end-50,z,s], cmap="viridis", clim=(0.1,cmax))
    n = s*2
    # ax.set_title("$n cells")
    ax.set_axis_off()

    # brain outline
    # @show size(im), size(small_outlines[50:end-100,1:end-50])
    outlines_to_show = small_outlines[yrange,xrange,z]
    # outlines_to_show = small_outlines[:,:,z]
    outlines_to_show = Float64.(outlines_to_show)
    if show_outline
        ax.imshow(outlines_to_show, alpha=outlines_to_show,cmap="binary_r")
    end
end

# axs[1].imshow(stim_masks[:, :,z,stimNum], cmap="gray")

# cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
# plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
#     norm=cnorm)
# may need to adjust if colorbar is cutoff
plt.tight_layout()
fig.subplots_adjust(right=0.93)
fig.subplots_adjust(wspace=0.04, hspace=0.01)
cbar_ax = fig.add_axes([0.94, 0.15, 0.0075, 0.7])
# cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
cbar = fig.colorbar(cim, cax=cbar_ax)
h = cbar.ax.set_ylabel("Δf/f",rotation=270,fontsize=7)
# h.set_rotation = 180
for (i,ax) in enumerate(axs)
    if i == 1
        ax.text(5,40,"1 target bilateral",fontsize=7, color="white")
    else
        ax.text(5,40,"$i",fontsize=7, color="white")
    end
end

um_per_px1 = tseries_units[2] / 1u"μm"
scalebar1 = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="lower right", color = "white",
    font_properties=Dict("size" => 7))
axs[1].add_artist(scalebar1)

# plotpath = joinpath(plot_dir,"$(analysis_name)_rdbu_titration_oneplane_df_f")
if no_df & show_outline
    plotpath = joinpath(plot_dir,"$(analysis_name)_outlined_titration_oneplane_df_f")
elseif show_outline
    plotpath = joinpath(plot_dir,"$(analysis_name)_viridis_medianfilt_outlined_titration_oneplane_df_f")
else
    plotpath = joinpath(plot_dir,"$(analysis_name)_viridis_medianfilt_titration_oneplane_df_f")
end
fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=300)
@show plotpath*".svg"
fig

#################################
## zoom-in on targeted cells
#################################
border = 5
stim_z = unique(stim_cells.z)
@assert length(stim_z) == 1
stim_z = stim_z[1]
stim_z = 5
# xrange = (minimum(stim_cells.x)-border):(maximum(stim_cells.x)+border)
# xrange = xrange .+ 5
xrange = 360:425
# yrange = (minimum(stim_cells.y)-border):(maximum(stim_cells.y)+border)
# yrange = Int(floor(yrange[1] + (yrange[end] - yrange[1])/2) + 15) : yrange[end]
yrange = 150:300
# yrange = 241:303

show_targets = false

N = 16
fig, axs = plt.subplots(2,8,figsize=(6,3))
axs = permutedims(axs,[2,1])
for i = 1:N
    ax = axs[i]
    stim_num = i
    stim_cells = unique(cells[cells.stimNum .== stim_num,:], [:x,:y])
    im = df_f[:,:,stim_z,stim_num]
    im = L.unsalt_n_pepper(df_f[:,:,stim_z,stim_num])

    ax.imshow(im[yrange,xrange],clim=(0,1.0))
    for (x,y,targetZ) in eachrow(stim_cells)
        cir_x = x - xrange[1]
        cir_y = y - yrange[1]
        circle = matplotlib.patches.Circle((cir_x,cir_y), targetSizePx/2,
            color="red",
            fill=false, lw=0.4)
        ax.add_patch(circle)
    end
    ax.set_axis_off()
end
plt.tight_layout()

plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_df_f_hab_target_closeup")
fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=300)
@show plotpath*".svg"

fig



# target_groups
##
tseries_subset = tseries[:,:,:,1:20:tseriesT];
maxTproj = mapslices(x->quantile(x[:],0.90), tseries_subset,dims=4)[:,:,:,1];
shifted_maxTproj = shift_every_other_row(maxTproj,-3);
##
fig, ax = plt.subplots()
z_im = zseries[yrange .* 2,xrange .* 2,imaging2zseries_plane[5]];
t_im = maxTproj[yrange,xrange,5];
t_im_shift = shifted_maxTproj[yrange,xrange,5];
# z_im = L.unsalt_n_pepper(z_im);
# z_im = adjust_histogram(imadjustintensity(z_im), GammaCorrection(0.9))
ax.imshow(z_im,cmap="Greys_r")
# ax.imshow(t_im,cmap="Greys_r")
# ax.imshow(t_im_shift,cmap="Greys_r")
for (x,y,targetZ) in eachrow(stim_cells)
    cir_x = x - xrange[1]
    cir_y = y - yrange[1]
    circle = matplotlib.patches.Circle((cir_x,cir_y), targetSizePx/2,
        color="red",
        fill=false, lw=0.4)
    ax.add_patch(circle)
end
ax.set_axis_off()

# plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_shifted_gray_hab_target_closeup")
plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_zseries_gray_hab_target_closeup")
fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=300)
@show plotpath
fig
##
save(joinpath(fish_dir, "zseries_imaging_planes.tif"), zseries[:,:,imaging2zseries_plane])
save(joinpath(fish_dir, "zseries.tif"), zseries)

##
dorsal_ventral = load(joinpath(fish_dir, "test_segs.nrrd"));
sum(dorsal_ventral[:,:,32])
Gray.(dorsal_ventral[:,:,32] .== 1) # dorsal
Gray.(dorsal_ventral[:,:,45] .== 2) # ventral