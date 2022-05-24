##
# ENV["DISPLAY"] = "localhost:12.0"
# ENV["DISPLAY"] = "/private/tmp/com.apple.launchd.5OQi0gJ6DL/org.xquartz:0"
##
# built from 2021-09-15_titration.jl. terrible hack for missing SLM trigger data
using ImageView
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation,
    Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
    Unitful, StatsBase
import PyPlot
import Plots
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
inset_axes = pyimport("mpl_toolkits.axes_grid1.inset_locator").inset_axes
matscale = pyimport("matplotlib_scalebar.scalebar")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman

matplotlib.rcParams["font.sans-serif"] = ["Arial", "sans-serif"]

# init_workers(16) # reduce memory usage
init_workers(8) # reduce memory usage
resources = Resources();
##
r = Recordings[
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
    # "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    # "2022-03-11_rschrmine_h2b6s_9dpf/fish1/TSeries-raphe-thalamus-118trial-002"
    # "2022-03-31_rschrmine_h2b6s_6dpf/fish1/TSeries-raphe-thalamus-118trial_no-autocorrect-350p-005"
    "2022-05-19_rschrmine_h2b6s_6dpf/fish1/TSeries-raphe-thalamus-5power-226trial-003"
](;
    resources...
    # LSTM tyh5 looks best
    # tseries_read_strategy = :lazy_hwzt,
    # tyh5_path="/data/dlab/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5",
    # tseries_dset=nothing
);

##
@pun tseries = r;
##
@pun (stim_start_idx, tseries) = r;
##

# snippet to recover from missing or bad VoltageRecording
green_bot = mean(tseries, dims=(1, 2, 3))[1, 1, 1, :];
# red_bot = loadBOT(tseries_dir, "Ch2");
##
Plots.plot(green_bot)
##
thresh = 6*std(green_bot)
putative_stim_starts = findall(green_bot .> thresh)
# not neighboring frames
valid_idx = vcat([1], diff(putative_stim_starts) .> 10) .== 1
putative_stim_starts = putative_stim_starts[valid_idx]
n_frames_between_stim = mode(diff(putative_stim_starts))
putative_stim_starts # no stim artifact for our control.... oops

##
# impute missing stimuli...
@pun trial_order = r
stim_start_idx = [putative_stim_starts[1]]
for i in 2:length(putative_stim_starts)
    isi = putative_stim_starts[i] - putative_stim_starts[i-1]
    prev = putative_stim_starts[i-1]
    # fill in missing stimuli
    while isi > 1.2 * n_frames_between_stim
        next = prev + n_frames_between_stim
        push!(stim_start_idx, next)
        isi -= n_frames_between_stim
        prev = next
    end
    push!(stim_start_idx, putative_stim_starts[i])
end

# missing first control stim & last control stim...
if length(stim_start_idx) == 116 && trial_order[1] == 3 && trial_order[end] == 3
    # missing first
    pushfirst!(stim_start_idx, putative_stim_starts[1]-n_frames_between_stim)
    # missing last
    push!(stim_start_idx, putative_stim_starts[end]+n_frames_between_stim)
end

##
# debug voltage recording
# @pun (voltageFile, tseriesZ) = r
# stimStartFrameIdx, stimEndFrameIdx, frameStartIdx, nstim_pulses, ttlStarts =
#     getStimTimesFromVoltages(voltageFile, tseriesZ;
#         stim_key="uncagingX") # 
# # stim_key="ell14_trigger") # 1
# # stim_key="AO0") # 0
# # stim_key="opto_gate") # 1
# # stim_key="OptogeneticPockels") # 78
# nstim_pulses
# ##
# import CSV
# voltages = CSV.File(open(read, voltageFile)) |> DataFrame
# rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])
# ##
# Plots.plot(voltages[10000:12000,
#     # "frame starts"]) # looks good
#     # "opto_gate"])
#     # "uncagingX"])
#     # "uncagingY"])
#     # "ell14_trigger"])
#     # "AO0"])
#     # "OptogeneticPockels"])
#     "FieldStimulator"])
##
@pun (tseriesH, tseriesW, tseriesZ, trial_order, vol_rate) = r
##
stim_end_idx = stim_start_idx .+ Int(ceil(vol_rate * 1))
pre = 10
post = 10
window_len = 10
trial_average = L.calc_trial_average(tseries, stim_start_idx,
    stim_end_idx, tseriesH, tseriesW, tseriesZ, trial_order;
    pre=pre, post=post)

##
# size(trial_average)
imshow(trial_average)
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
##
@pun (tseries_read_strategy, exp_name, fish_dir) = r
if occursin("tyh5", string(tseries_read_strategy))
    avgstim_path = joinpath(fish_dir, exp_name * "_lstm_avgStim.h5")
elseif occursin("tiff", string(tseries_read_strategy))
    avgstim_path = joinpath(fish_dir, exp_name * "_raw_avgStim.h5")
elseif occursin("hwzt", string(tseries_read_strategy))
    avgstim_path = joinpath(fish_dir, exp_name * "_kalman_avgStim.h5")
end


try
    if occursin("tyh5", string(tseries_read_strategy))
        global analysis_name = "lstm_imap"
        global trial_average = h5read(avgstim_path, "/lstm")
    elseif occursin("hwzt", string(tseries_read_strategy))
        global analysis_name = "kalman_imap"
        global trial_average = h5read(avgstim_path, "/kalman")
    else
        global analysis_name = "raw_imap"
        # global trial_average = h5read(avgstim_path, "/lstm");
        global trial_average = h5read(avgstim_path, "/raw")
    end
catch
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println(stdout)
    end
    println("couldn't read $avgstim_path")
    global trial_average = r[:trial_average]
    if tseries_read_strategy == :lazy_tyh5
        h5write(avgstim_path, "/lstm", trial_average)
    elseif tseries_read_strategy == :lazy_tiff
        h5write(avgstim_path, "/raw", trial_average)
    elseif tseries_read_strategy == :lazy_hwzt
        h5write(avgstim_path, "/kalman", trial_average)
    end
end

size(trial_average)
##
imap = influence_map(trial_average, window_len);
sum(imap, dims=[1, 2, 3])[1, 1, 1, :]
##
# imshow(imap)
##

s = stim_start_idx[1]
# imshow(shift_every_other_row(tseries[:,:,:,s-10:s+15],-3))
##
@pun (exp_date, lateral_unit, nStimuli, etl_vals, recording_folder, fish_name, plot_dir) = r
etl_vals = collect(-40:10:50)
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
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)
global fig = plt.figure()
@warn "df_f denominator ϵ may have changed"
for stimNum in 1:nStimuli
    # for stimNum in 16:nStimuli
    df_f = imap[:, :, :, stimNum]
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    if tseriesZ > 5
        figW, figH = (figB * tseriesZ / 2, figB * 2)
        global fig = plt.figure(figsize=(figW, figH))
        fig, axs = plt.subplots(2, Int(tseriesZ / 2), figsize=(figW, figH))
    else
        figW, figH = (figBtseriesZZ, figB)
        global fig = plt.figure(figsize=(figW, figH))
        fig, axs = plt.subplots(1, Z, figsize=(figW, figH))
    end
    if tseriesZ == 1
        axs = [axs]
    end
    # ax = axs[2]
    # ax = axs
    for z in 1:tseriesZ
        ax = axs[z]
        global cim = ax.imshow(df_f[:, :, z], cmap="RdBu_r",
            norm=cnorm)
        ax.set_axis_off()
        ax.set_title("$(Int(round(etl_vals[z],digits=0)))μm")
        # this will make extra circles (1 extra per repetition...)
        # for (x, y, targetZ) in eachrow(unique(cells[cells.stimNum.==stimNum, [:x, :y, :z]]))
        #     if z == targetZ
        #         circle = matplotlib.patches.Circle((x, y), targetSizePx, color="k",
        #             fill=false, lw=0.4, alpha=0.3)
        #         ax.add_patch(circle)
        #     end
        # end
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
    path = joinpath(plot_dir, "$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_stim$stimNum")
    @show path * ".svg"
    fig.savefig(path * ".svg", dpi=600)
    fig.savefig(path * ".png", dpi=600)
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

    # to manually draw:
    # plot once with uncommented, trace diencephalon in illustrator, re-run with commented
    # "Diencephalon -"
]
region_masks = [L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane, region; outline=false) for region in REGION_LIST];
# dilate on outside for thicker line

##
region_outlines = [dilate(L.mask2outline(L.unsalt_n_pepper(regm[2]; felz_min=5000)), [1, 2])
                   for regm in region_masks];
# for regm in region_masks[[1]]];
outlines = cat(region_outlines..., dims=4);
outlines = maximum(outlines, dims=4)[:, :, :, 1];
small_outlines = imresize(outlines, (tseriesH, tseriesW, tseriesZ));
Gray.(outlines[:, :, 3])

#################################
## Titration one-plane
#################################
figB = 1.6

nStimuli = size(trial_average)[end-1]
# nStimuli = size(long_avg)[end-1]
volRate = 3
max_frames = size(trial_average, 5)
Z = size(trial_average, 3)
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)

cmax = 1.0
# cmax = 2
# cmin = -1.0
cmin = 0
# cmin = -0.5
nseconds = 5
pre = Int(ceil(nseconds * volRate)) + 1
post = Int(ceil(nseconds * volRate)) + 1

window = minimum([Int(ceil(3 * volRate)), max_frames])
@assert (window <= post) & (window <= pre)
# @assert Z == 1

# cmax = 4
# cmin = -0.75
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
stimNum = 16

# f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
# TODO: should window be immediately after stim..? right now it's
# ~2 seconds after stim...?
f = mean(trial_average[:, :, :, :, end-window+1:end], dims=5)[:, :, :, :, 1]
f0 = mean(trial_average[:, :, :, :, 1:window], dims=5)[:, :, :, :, 1]
df = f - f0
df_f = df ./ f0
nplots = 15
@assert nplots == nStimuli - 1

# figW,figH = (figB*4, figB*4)
# two-column width (mm)
figB = Float64(uconvert(u"inch", (183 - 3)u"mm") / 4u"inch")
figW, figH = (figB * 4, figB * 3)
##

fig = plt.figure(figsize=(figW, figH),
    constrained_layout=false)
# constrained_layout=true)
# each titration panel is 6x6
gridbase = 6
first_col_width = 4
gs = fig.add_gridspec(gridbase * 4, gridbase * 4 + first_col_width)
# ax1 = fig.add_subplot(gs[1:6, 1]) # 1-cell stim close-up
ax1 = fig.add_subplot(py"$(gs)[0:2*$gridbase,0:$(first_col_width)]") # 1-cell stim close-up
# ax2 = fig.add_subplot(gs[7:12, 1]) # colorbar
ax2 = fig.add_subplot(py"$(gs)[2*$gridbase:,0:1]") # colorbar
axs = [ax1, ax2]
for (i, j) in Iterators.product(1:4, 1:4)
    # first at `gridbase`, second at `2*gridbase`...
    is = (i - 1) * gridbase + first_col_width
    ie = i * gridbase + first_col_width
    # first at 1, second at 4
    js = (j - 1) * gridbase
    je = j * gridbase
    # println(py"($js,$je,$is,$ie)")
    push!(axs, fig.add_subplot(py"$(gs)[$js:$je,$is:$ie]"))
end

# function example_plot(ax, fontsize=12)
#     ax.plot([1, 2])

#     ax.locator_params(nbins=3)
#     ax.set_xticklabels([])
#     ax.set_yticklabels([])
# end
# for ax in axs
#     example_plot(ax)
# end
# fig
z = 3
felz_k = 50
felz_min = 10
# show_outline = false
no_df = false
show_outline = true
# no_df = true
# yrange = 65:tseriesH-100
yrange = 85:tseriesH-120
xrange = 70:tseriesW-90
zoom_xrange = 360:425
zoom_yrange = (150:300) .+ 10
inset1_xrange = (1:40) .+ 133
inset1_yrange = (256:296) .- 75
inset2_xrange = (1:40) .+ 133
inset2_yrange = (256:296) .- 5
cmap = "inferno"
# cmap = "magma"
# cmap = "viridis"
lw = 0.5
function style_inset(ax, color="white")
    ax.spines["bottom"].set_color(color)
    ax.spines["top"].set_color(color)
    ax.spines["right"].set_color(color)
    ax.spines["left"].set_color(color)
    ax.set_xticks([])
    ax.set_xticks([], minor=true)
    ax.set_yticks([])
    ax.set_yticks([], minor=true)
end

um_per_px1 = tseries_units[2] / 1u"μm"

# xrange = 40:tseriesW-70
for (s, ax) in zip(1:nStimuli, axs[3:end])
    im = df_f[:, :, z, s]
    im = opening_median(im)
    segments = felzenszwalb(im, felz_k, felz_min)
    im[segments.image_indexmap.==1] .= 0
    main_im = im[yrange, xrange]
    # global cim = ax.imshow(df_f[50:end-100,1:end-50,z,s], cmap="RdBu_r",
    # global cim = ax.imshow(im, cmap="RdBu_r",
    if ~no_df
        global cim = ax.imshow(main_im, cmap=cmap,
            clim=(cmin, cmax), interpolation="none")
        # norm=cnorm, interpolation="none")
    end
    # ax.imshow(outlines[:,:,z], alpha=outlines[:,:,z])
    # global cim = ax.imshow(df_f[50:end-100,1:end-50,z,s], cmap=cmap, clim=(0.1,cmax))
    n = s * 2
    # ax.set_title("$n cells")
    ax.set_axis_off()

    # brain outline
    # @show size(im), size(small_outlines[50:end-100,1:end-50])
    outlines_to_show = small_outlines[yrange, xrange, z]
    # outlines_to_show = small_outlines[:,:,z]
    outlines_to_show = Float64.(outlines_to_show)
    if show_outline
        if s == 1
            outlines_to_show[1:100, 1:30] .= 0
        end
        ax.imshow(outlines_to_show, alpha=outlines_to_show, cmap="binary_r")
    end

    if s == 1
        # add close-up
        axs[1].imshow(im[zoom_yrange, zoom_xrange], clim=(0, 1.0), cmap=cmap)
        axs[1].set_axis_off()
        if show_outline
            outlines_to_show = small_outlines[zoom_yrange, zoom_xrange, z]
            outlines_to_show = Float64.(outlines_to_show)
            # axs[1].imshow(outlines_to_show, alpha=outlines_to_show,cmap="binary_r")
        end
        rect = matplotlib.patches.Rectangle(
            (1, 1),
            # length(zoom_xrange)-3, length(zoom_yrange)-3,
            length(zoom_xrange), length(zoom_yrange) - 3,
            linestyle="--", linewidth=lw, edgecolor="w", facecolor="none")
        axs[1].add_patch(rect)


        # dotted rect for big fig
        rect = matplotlib.patches.Rectangle(
            (zoom_xrange[1] - xrange[1], zoom_yrange[1] - yrange[1]),
            zoom_xrange[end] - zoom_xrange[1], zoom_yrange[end] - zoom_yrange[1],
            # zoom_xrange[end]-zoom_xrange[1]-5, zoom_yrange[end]-zoom_yrange[1],
            linestyle="--", linewidth=lw, edgecolor="w", facecolor="none")
        ax.add_patch(rect)
    else
        # add inset
        # upper right
        rect1 = matplotlib.patches.Rectangle(
            (inset1_xrange[1] - xrange[1], inset1_yrange[1] - yrange[1]),
            length(inset1_xrange), length(inset1_yrange),
            linestyle="--", linewidth=lw, edgecolor="w", facecolor="none")
        ax.add_patch(rect1)

        rect2 = matplotlib.patches.Rectangle(
            (inset2_xrange[1] - xrange[1], inset2_yrange[1] - yrange[1]),
            length(inset2_xrange), length(inset2_yrange),
            linestyle="--", linewidth=lw, edgecolor="w", facecolor="none")
        ax.add_patch(rect2)


        axins1 = inset_axes(ax, height="100%", width="100%",
            bbox_to_anchor=(0.13, 0.76, 0.28, 0.28), bbox_transform=ax.transAxes, loc="upper right")

        axins1.imshow(im[inset1_yrange, inset1_xrange], clim=(0, 1.0), cmap=cmap)

        style_inset(axins1)
        # lower right
        axins2 = inset_axes(ax, height="100%", width="100%",
            bbox_to_anchor=(0.13, -0.04, 0.28, 0.28), bbox_transform=ax.transAxes, loc="lower right")

        # make axis a dotted white
        # for spine in ["right", "top", "left", "bottom"]
        #     axins1.spines[spine].set_linestyle((0,(4,3)))
        #     axins2.spines[spine].set_linestyle((0,(4,3)))
        # end


        axins2.imshow(im[inset2_yrange, inset2_xrange], clim=(0, 1.0), cmap=cmap)
        style_inset(axins2)
        if s == 2
            # scalebar for second panel only
            scalebar = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.25, box_alpha=0,
                scale_loc="left", location="lower right", color="white",
                pad=0.1, font_properties=Dict("size" => 7),
                scale_formatter=py"""lambda value, unit: "" """)
            axins2.add_artist(scalebar)
        end

        # top
        ls = (0, (3, 2))
        x1 = inset1_xrange[1] - xrange[1]
        y1 = inset1_yrange[1] - yrange[1]
        x2 = inset1_xrange[end] - xrange[1]
        y2 = inset1_yrange[end] - yrange[1]
        ax.plot([x1, x1 - 25], [y1, y1 - 9], ls=ls, color="white", lw=0.5)
        ax.plot([x2, x2 + 24], [y1, y1 - 9], ls=ls, color="white", lw=0.5)

        # bot
        x1 = inset2_xrange[1] - xrange[1]
        y1 = inset2_yrange[1] - yrange[1]
        x2 = inset2_xrange[end] - xrange[1]
        y2 = inset2_yrange[end] - yrange[1]
        ax.plot([x1, x1 - 25], [y2, y2 + 13], ls=ls, color="white", lw=0.5)
        ax.plot([x2, x2 + 22], [y2, y2 + 13], ls=ls, color="white", lw=0.5)

        # annotations
        if s >= 14
            # top inset
            circ1 = matplotlib.patches.Circle(
                (32, 32), 5,
                linestyle="-", linewidth=lw, edgecolor="w", facecolor="none")
            axins1.add_patch(circ1)
        end

        if s >= 3
            # bottom inset
            circ2 = matplotlib.patches.Circle(
                (17, 21), 5,
                linestyle="-", linewidth=lw, edgecolor="w", facecolor="none")
            axins2.add_patch(circ2)
        end

    end
end

# axs[1].imshow(stim_masks[:, :,z,stimNum], cmap="gray")

# cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
# plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
#     norm=cnorm)
# may need to adjust if colorbar is cutoff
plt.tight_layout()
# fig.subplots_adjust(right=0.93)
fig.subplots_adjust(wspace=0.04, hspace=0.04)
# cbar_ax = fig.add_axes([0.94, 0.15, 0.0075, 0.7])
cbar = fig.colorbar(cim, cax=axs[2], aspect=20)
h = cbar.ax.set_ylabel("Δf/f", rotation=270, fontsize=7)
for (i, ax) in enumerate(axs[3:end])
    if i == 1
        ax.text(3, 80, "1\ntarget\nbilateral", fontsize=7, color="white")
    else
        ax.text(3, 25, "$i", fontsize=7, color="white")
    end
end


scalebars = [matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="lower right", color="white",
    font_properties=Dict("size" => 7),
    scale_formatter=py"""lambda value, unit: "" """)
             for _ in 1:2]
axs[3].add_artist(scalebars[1])
axs[1].add_artist(scalebars[2])

# plotpath = joinpath(plot_dir,"$(analysis_name)_rdbu_titration_oneplane_df_f")
if no_df & show_outline
    plotpath = joinpath(plot_dir, "$(analysis_name)_outlined_titration_oneplane_df_f")
elseif show_outline
    plotpath = joinpath(plot_dir, "$(analysis_name)_$(cmap)_medianfilt_outlined_titration_oneplane_df_f")
else
    plotpath = joinpath(plot_dir, "$(analysis_name)_$(cmap)_medianfilt_titration_oneplane_df_f")
end
fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir, "$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir, "$(plotpath).pdf"),
    dpi=300)
@show plotpath * ".svg"
fig
