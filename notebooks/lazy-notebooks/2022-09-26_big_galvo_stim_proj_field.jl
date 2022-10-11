##
ENV["DISPLAY"] = "localhost:10.0"
# ENV["DISPLAY"] = "/private/tmp/com.apple.launchd.5OQi0gJ6DL/org.xquartz:0"
##
# use for missing voltage / trigger data for using alternate galvo stim path
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation,
    Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
    Unitful, StatsBase, LibExpat, ImageDraw, Arrow
using ImageView
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

init_workers(16) # reduce memory usage
# init_workers(8) # reduce memory usage
resources = Resources();
##
r = Recordings[
# no SLM trigger since using laser stim
# "2022-09-14_rschrmine_h2b6s_6dpf/fish2_chrmine-neg/TSeries_stim-lthal_10z_3-32-aperature-055"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish2_chrmine-neg/TSeries_stim-lthal_10z_1-aperature-053"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish2_chrmine-neg/TSeries_stim-lthal_10z_3-32-aperature-052"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish2_chrmine-neg/TSeries_stim-lthal_10z_1-aperature-050"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_3-32-aperature_4MHz-077"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_1-aperature_4MHz-076"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_3-32-aperature_4MHz-074"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_1-aperature_4MHz-073"
# "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_3-32-aperature_4MHz-072"
    "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_1-aperature_4MHz-070"
](;
    resources...,
    tseries_read_strategy=:tiff
);
flip_x = occursin(r"fish2_chrmine-neg", r[:uri]) # only flip fish2
analysis_name = occursin(r"1-aperature", r[:uri]) ? "big-NA" : "small-NA"
fwhm_val = analysis_name == "big-NA" ? 1.3 : 7.3
#
@pun (tseries, tseries_dir, tseries_units) = r;
markpoints_xml = glob_one_file("*_MarkPoints.xml", tseries_dir);
big_target_center = read_markpoints_series(markpoints_xml)[1]
mp_xml = parseXML(markpoints_xml)
spiral_μm = 1.0u"μm" * parse(Int, mp_xml[xpath"/PVMarkPointSeriesElements/PVMarkPointElement/PVGalvoPointElement/Point/@SpiralSizeInMicrons"][1])
spiral_px = spiral_μm / tseries_units[1]


# snippet to recover from missing or bad VoltageRecording
green_bot = mean(tseries, dims=(1, 2, 3))[1, 1, 1, :];
# red_bot = loadBOT(tseries_dir, "Ch2");
# Plots.plot(green_bot)

thresh = 5 * std(green_bot)
putative_stim_starts = findall(green_bot .> thresh)
# not neighboring frames
valid_idx = vcat([1], diff(putative_stim_starts) .> 10) .== 1
putative_stim_starts = putative_stim_starts[valid_idx]
n_frames_between_stim = mode(diff(putative_stim_starts))
putative_stim_starts # no stim artifact for our control.... oops


# impute missing stimuli...
stim_start_idx = [putative_stim_starts[1]]
for i in 2:length(putative_stim_starts)
    isi = putative_stim_starts[i] - putative_stim_starts[i-1]
    prev = putative_stim_starts[i-1]
    # fill in missing stimuli
    while isi > 1.8 * n_frames_between_stim
        next = prev + n_frames_between_stim
        push!(stim_start_idx, next)
        isi -= n_frames_between_stim
        prev = next
    end
    push!(stim_start_idx, putative_stim_starts[i])
end

@assert abs(maximum(diff(stim_start_idx)) - minimum(diff(stim_start_idx))) < 3

@pun (tseriesH, tseriesW, tseriesZ, vol_rate) = r
nTrials = length(stim_start_idx)
stim_end_idx = stim_start_idx .+ Int(ceil(vol_rate * 1)) .+ 1
pre = 10
post = 10
window_len = 10
trial_average = L.calc_trial_average(tseries, stim_start_idx,
    stim_end_idx .+ 1, tseriesH, tseriesW, tseriesZ, ones(Int, length(stim_start_idx));
    pre=pre, post=post);
# imshow(trial_average[:, :, :, 1, :])

imap = influence_map(trial_average, window_len - 1);


@pun (etl_vals, fish_name, exp_name, recording_folder, plot_dir) = r
nStimuli = 1
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)
window = window_len - 1
with_circles = true

cmax = 2.0
cmin = -0.2

ϵ = 0.0
# ϵ = 0.05

cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)
global fig = plt.figure()
@warn "df_f denominator ϵ may have changed"
for stimNum in 1:nStimuli
    df_f = imap[:, :, :, stimNum]
    df_f = L.unsalt_n_pepper(df_f)
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    if tseriesZ % 2 == 1
        tseriesZ_display = tseriesZ - 1
    else
        tseriesZ_display = tseriesZ
    end
    if tseriesZ > 5
        figW, figH = (figB * tseriesZ_display / 2, figB * 2)
        global fig = plt.figure(figsize=(figW, figH))
        fig, axs = plt.subplots(2, Int(tseriesZ_display / 2), figsize=(figW, figH))
    else
        figW, figH = (figBtseriesZZ, figB)
        global fig = plt.figure(figsize=(figW, figH))
        fig, axs = plt.subplots(1, Z, figsize=(figW, figH))
    end
    if tseriesZ_display == 1
        axs = [axs]
    end
    # ax = axs[2]
    # ax = axs
    for z in 1:tseriesZ_display
        ax = axs[z]
        global cim = ax.imshow(df_f[:, :, z], cmap="RdBu_r",
            norm=cnorm)
        ax.set_axis_off()
        z_um = Int(round(etl_vals[z], digits=0))
        ax.set_title("$(z_um)μm")


        alpha = z_um == 0 ? 0.7 : 0.2
        circle = matplotlib.patches.Circle((big_target_center[2], big_target_center[1]), spiral_px / 2, color="k",
            fill=false, lw=0.4, alpha=alpha)
        ax.add_patch(circle)
        if flip_x
            ax.invert_xaxis()
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
    path = joinpath(plot_dir, "$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_stim$stimNum")
    if with_circles
        path *= "_circle-targets"
    end
    @show path * ".svg"
    fig.savefig(path * ".svg", dpi=600)
    fig.savefig(path * ".png", dpi=600)
    fig
end
fig

# save arrow file for ROI by Z

mask = Gray.(zeros(Bool, tseriesH, tseriesW))
for target in [big_target_center]
    targetSize = spiral_px / 2
    draw!(mask, Ellipse(CirclePointRadius(target[2], target[1], targetSize)))
end
mask
mask = findall(channelview(mask))

#
ϵ = 1.0
df = DataFrame()
nskip = 2
# assumes only one stimulation type
for t in 1:nTrials
    sidx = stim_start_idx[t]
    eidx = stim_end_idx[t]
    for z in 1:tseriesZ
        start = sidx - pre
        stop = eidx + post
        f0 = tseries[:, :, z, start:sidx-nskip]
        f = tseries[:, :, z, eidx+nskip:stop]
        f0 = mean(f0[mask, :])
        f = mean(f[mask, :])
        df_f = (f - f0) / (f0 + ϵ)
        df = vcat(df,
            DataFrame(
                "Δf/f" => df_f,
                "trial" => t,
                "etl" => etl_vals[z],
                "recording_folder" => recording_folder,
                "fish_name" => fish_name,
                "exp_name" => exp_name,
                "fwhm" => fwhm_val
            )
        )
    end
end
df


arrow_df_path = tseries_dir * "_roiDF.arrow"
open(arrow_df_path, "w") do io
    Arrow.write(io, df)
end
println("saved arrow file: $arrow_df_path")

# plot one file...
rdf = Arrow.Table(arrow_df_path) |> DataFrame


points = Gadfly.layer(rdf, x="etl", y="Δf/f", color="fwhm", Geom.point)
zmean = combine(groupby(rdf, "etl"), "Δf/f" => mean, "fwhm" => mean, renamecols=false)
line = Gadfly.layer(zmean, x="etl", y="Δf/f", color="fwhm", Geom.line)
Gadfly.plot(points, line)


## plot all files
root_dir = "/data/dlab/b118/"
chrmine_df_paths = glob("*.arrow", joinpath(root_dir, "2022-09-14_rschrmine_h2b6s_6dpf/fish3/"))
control_df_paths = glob("*.arrow", joinpath(root_dir, "2022-09-14_rschrmine_h2b6s_6dpf/fish2_chrmine-neg/"))

chr_df = DataFrame()
for p in chrmine_df_paths
    rdf = Arrow.Table(p) |> DataFrame
    chr_df = vcat(chr_df, rdf)
end
chr_df

con_df = DataFrame()
for p in control_df_paths
    rdf = Arrow.Table(p) |> DataFrame
    con_df = vcat(con_df, rdf)
end
con_df

##
points = Gadfly.layer(chr_df, x="etl", y="Δf/f", color="fwhm", Geom.beeswarm)
zmean = combine(groupby(chr_df, [:etl, :fwhm]), "Δf/f" => mean, renamecols=false)
line = Gadfly.layer(zmean, x="etl", y="Δf/f", color="fwhm", Geom.line)
Gadfly.plot(points, line)

##
p1 = Gadfly.plot(chr_df, xgroup="fwhm", Geom.subplot_grid(
    Gadfly.layer(chr_df, x="etl", y="Δf/f", color="fwhm", Geom.beeswarm),
    Gadfly.layer(zmean, x="etl", y="Δf/f", color="fwhm", Geom.line)
))

##
p1 = Gadfly.plot(chr_df, xgroup="fwhm", x="etl", y="Δf/f",
    Geom.subplot_grid(Geom.boxplot))
##
p1 = Gadfly.plot(chr_df, x="etl", y="Δf/f", color=:fwhm,
    Scale.x_discrete(levels=etl_vals),
    Guide.yticks(ticks=[-0.25, 0., 0.25, 0.5, 0.75]),
    Coord.cartesian(ymin=-0.3, ymax=0.8),
    Scale.color_discrete(),
    Guide.xlabel("etl (μm)"),
    Geom.boxplot(suppress_outliers=true))

p2 = Gadfly.plot(con_df, x="etl", y="Δf/f", color=:fwhm,
    Scale.x_discrete(levels=etl_vals),
    Guide.yticks(ticks=[-0.25, 0., 0.25, 0.5, 0.75]),
    Coord.cartesian(ymin=-0.3, ymax=0.8),
    Scale.color_discrete(),
    Guide.xlabel("etl (μm)"),
    Geom.boxplot(suppress_outliers=true))

# p1_path = joinpath(splitpath(chrmine_df_paths[1])[1:end-1]..., "raw/plots/na_comparison.svg")
# img1 = SVG(p1_path, 6inch, 5inch)
# p2_path = joinpath(splitpath(control_df_paths[1])[1:end-1]..., "raw/plots/na_comparison.svg")
# img2 = SVG(p2_path, 6inch, 5inch)

p1_path = joinpath(splitpath(chrmine_df_paths[1])[1:end-1]..., "raw/plots/na_comparison.png")
img1 = PNG(p1_path, 6inch, 5inch)
p2_path = joinpath(splitpath(control_df_paths[1])[1:end-1]..., "raw/plots/na_comparison.png")
img2 = PNG(p2_path, 6inch, 5inch)

Gadfly.draw(img1, p1)
Gadfly.draw(img2, p2)

p1_path, p2_path
##
p1
# Gadfly.vstack(p1, p2)