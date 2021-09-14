# ENV["DISPLAY"] = "localhost:11.0"
##
# using ImageView
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Unitful, ImageDraw
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
# init_workers()
##
resources = Resources();
recording1 = Recordings[
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-lrhab-118trial-060"
    "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061"
](;resources...,
    tseries_read_strategy = :lazy_hwzt,
    tseries_suffix = "kalman",
    tseries_dset = nothing,
    tyh5_path = "/data/dlab/b115/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061_kalman.h5"
);
# second-order experiment
recording2 = Recordings[
    # actually from 6-02...
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
    "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
](;resources...,
    tseries_read_strategy = :lazy_hwzt,
    tseries_suffix = "kalman",
    tseries_dset = nothing,
    tyh5_path = "/data/dlab/b115/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072_kalman.h5"
);
##
@pun zbrain_masks = resources;
@pun (trial_average,vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, stim_end_idx, stim_start_idx, nstim_pulses, window_secs,
    fish_dir, tseries_units, lateral_unit, tseries_read_strategy, trial_average_path,
    slm_num, zOffset, tseries, tseriesT, tseriesH, tseriesW
) = recording2;
##
@pun (zseries, imaging2zseries_plane, zseries_xml, zseries_zaxes,
    zseries_dir) = recording2;
##
@show nstim_pulses
@assert nstim_pulses == 10
tseries_read_strategy
##
hab_trial_average = recording1[:trial_average];
hab_trial_average_dset = recording1[:trial_average_dset];
hab_trial_average_path = recording1[:trial_average_path];
hab_cells = recording1[:cells];
hab_nStimuli = recording1[:nStimuli];
hab_target_groups = recording1[:target_groups];
##
hab_imap = influence_map(hab_trial_average, window_len);
imap = influence_map(trial_average, window_len);
## draw regions example

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

# cmax = 3
cmax = 1
cmin = 0
# cmin = -0.2

ϵ = 0.0
# ϵ = 0.05

# cmax = 0.5
# cmin = -0.1
# cmax = 4
# cmin = -0.75
# cmin = -4
ncol = 3
nrows = 5
# figB = Float64(uconvert(u"inch", 183u"mm") / (ncol*u"inch"))
# figW,figH = (figB*ncol, figB*3)
figH = Float64(uconvert(u"inch", 13u"cm") / u"inch")
figW = 9.3/10 * figH

fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(nrows, 3, figsize=(figW,figH))

ztotext(z) = "$(Int(round(etl_vals[z],digits=0)))"

"Remove salt & pepper."
function prettify_img(im)
    im = opening_median(im)
    segments = felzenszwalb(im, felz_k, felz_min)
    im[segments.image_indexmap .== 1] .= 0
    im
end
zplanes = [5,10]
yrange = 40:400
xrange = 1:512
felz_k = 50
felz_min = 10

n = 100
zoom_yrange = 200:200+n
zoom_xrange = 175:175+n


target_planes = unique(vcat(cells.z, hab_cells.z))
ntarget_planes = length(target_planes)
seeds = [RGB(1,1,1)]
colors = distinguishable_colors(ntarget_planes, seeds, dropseed=true)
plane_color = Dict(z => c for (z,c) in zip(target_planes, colors))



# make each stim a row...
# for row in 1:(nStimuli+hab_nStimuli-1)
for row in [3,1,2,4,5]
    if row >3
        df_f = imap[:,:,:,row-3]
        stimNum = row - 3
        stim_cells = cells
    else
        df_f = hab_imap[:,:,:,row]
        stim_cells = hab_cells
        stimNum = row
    end
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)

    # show targets on first column
    ax = axs[row]
    theSize = collect(size(df_f)[1:2])
    theSize[1] = length(yrange)
    ax.imshow(zeros(theSize...), cmap="RdBu_r", clim=(-1,1))
    each_target = unique(stim_cells[stim_cells.stimNum .== stimNum,[:x,:y,:z]])

    for (x,y,targetZ) in eachrow(each_target)
        circle = matplotlib.patches.Circle((x,y-yrange[1]), targetSizePx,
            color=rgb2tuple(plane_color[targetZ]),
            fill=false, lw=0.2)
        ax.add_patch(circle)
    end
    ax.set_axis_off()
    
    for (j,z) in enumerate(zplanes)
        ax_num = nrows * j + row
        ax = axs[ax_num]
        im = df_f[:,:,z]
        im = opening_median(im)
        segments = felzenszwalb(im, felz_k, felz_min)
        im[segments.image_indexmap .== 1] .= 0

        global cim = ax.imshow(im[yrange,:], cmap="viridis",
            clim=(cmin,cmax), interpolation="none")
        ax.text(5,yrange[end]-50,"z = $(Int(etl_vals[z]))μm",fontsize=7, color="white")
        # zoom in ROI (compile in seperate plot)
        if ((row==1) || (row==4)) && (j==1)
            ls = row == 1 ? "--" : "-"
            @show  (zoom_xrange[1] - xrange[1], zoom_yrange[1] - yrange[1]),
                zoom_xrange[end]-zoom_xrange[1], zoom_yrange[end]-zoom_yrange[1]
            rect = matplotlib.patches.Rectangle(
                (zoom_xrange[1] - xrange[1], zoom_yrange[1] - yrange[1]),
                zoom_xrange[end]-zoom_xrange[1], zoom_yrange[end]-zoom_yrange[1],
                linestyle=ls, linewidth=1, edgecolor="w", facecolor="none")
            ax.add_patch(rect)
        end

        ax.set_axis_off()
    end
end
sorted_plane_color = sort(collect(plane_color);lt=(a,b)->isless(a[1],b[1]))
patch_handles = [ matplotlib.lines.Line2D([0],[0],
    markerfacecolor=nothing, markeredgecolor=rgb2tuple(c), color="w",
    label="$(Int(etl_vals[z]))μm", marker="o", markersize=Int(round(targetSizePx))-2)
    for (z,c) in sorted_plane_color]

um_per_px1 = tseries_units[2] / 1u"μm"
scalebar1 = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="upper right", font_properties=Dict("size" => 7), color="w")
axs[6].add_artist(scalebar1)

plt.tight_layout()
fig.subplots_adjust(wspace=0.01, hspace=0.01)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.91, 0.09, 0.0075, 0.7])
# cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
cbar = fig.colorbar(cim, cax=cbar_ax)
h = cbar.ax.set_ylabel("Δf/f",rotation=270, fontsize=7)


legend = axs[1].legend(handles=patch_handles, title="target z-plane", prop=Dict("size" => 6),
    # loc="lower left", bbox_to_anchor=(0.0, 0.0))
    loc="lower left", bbox_to_anchor=(0.0, 0.0))
plt.setp(legend.get_title(),fontsize=6)


path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_8panel")
@show path*".svg"
fig.savefig(path*".svg", dpi=300)
fig.savefig(path*".png", dpi=300)
fig.savefig(path*".pdf", dpi=300)

fig
## get max t project
tseries_subset = tseries[:,:,:,1:20:tseriesT];
maxTproj = mapslices(x->quantile(x[:],0.90), tseries_subset,dims=4)[:,:,:,1];

## detail zoom
z = zplanes[1]
# z = 5
stimNum = 1
figH = Float64(uconvert(u"inch", 8u"cm") / u"inch")
figW = figH / 2

fig, axs = plt.subplots(2,1, figsize=(figW, figH))
# zoom_yrange = 1:512
# zoom_xrange = 1:512
im = hab_imap[:,:,z,1]
# im = mean(hab_trial_average,dims=[4,5])[:,:,z,1,1]
im = opening_median(im)
segments = felzenszwalb(im, felz_k, felz_min)
im[segments.image_indexmap .== 1] .= 0

global cim = axs[1].imshow(im[zoom_yrange,zoom_xrange], cmap="viridis",
    clim=(cmin,cmax), interpolation="none")

z_im = zseries[zoom_yrange*2,zoom_xrange*2,imaging2zseries_plane[z]]
# z_im = zseries[zoom_yrange*2,zoom_xrange*2,24]
proj_im = imap[:,:,z,stimNum]
proj_im = opening_median(proj_im)
segments = felzenszwalb(proj_im, felz_k, felz_min)
proj_im[segments.image_indexmap .== 1] .= 0

# axs[2].imshow(z_im)
axs[2].imshow(proj_im[zoom_yrange,zoom_xrange], cmap="viridis",
    clim=(cmin,cmax*2), interpolation="none")

struct_im = maxTproj[zoom_yrange,zoom_xrange,z]
# struct_im = adjust_histogram(struct_im, GammaCorrection(0.5));
# struct_im = adjust_histogram(struct_im, ContrastStretching(t=0.2,slope=1.5))
# struct_im = adjust_histogram(struct_im, ContrastStretching(t=0.5, slope=0.1))
local_win_px = 7u"μm"/lateral_unit
nblocks = Int(round(size(struct_im,1) / local_win_px))
eq = ImageContrastAdjustment.AdaptiveEqualization(nbins = 64, rblocks = nblocks, cblocks = nblocks, clip = 0.5)
struct_im = imadjustintensity(struct_im)
struct_im = adjust_histogram(struct_im, eq)

# axs[3].imshow(struct_im, cmap="gray",
#     interpolation="none")


plane_stim_cells = unique(cells[(cells.stimNum .== stimNum) .& (cells.z .== z),
    [:x,:y]])

for (x,y) in eachrow(plane_stim_cells)
    circle = matplotlib.patches.Circle((x-zoom_xrange[1],y-zoom_yrange[1]), targetSizePx,
        color="red",
        fill=false, lw=0.5)
        axs[1].add_patch(circle)
    circle = matplotlib.patches.Circle((x-zoom_xrange[1],y-zoom_yrange[1]), targetSizePx,
        color="red",
        fill=false, lw=0.5)
    axs[2].add_patch(circle)
end
axs[1].set_axis_off()    
axs[2].set_axis_off()    
# axs[3].set_axis_off()    

plt.tight_layout()
fig.subplots_adjust(wspace=0.01, hspace=0.01)

scalebar2 = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", color="white", location="lower left", font_properties=Dict("size" => 7))
axs[1].add_artist(scalebar2)



path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_zoomdetail")
@show path*".svg"
fig.savefig(path*".svg", dpi=300)
fig.savefig(path*".png", dpi=300)
fig.savefig(path*".pdf", dpi=300)

fig


## plot overlap

medulla_mask = Gray.(zeros(Bool,size(imap)[1:2]))
draw!(medulla_mask, Ellipse(
    Point(235,238), 60, 60))
medulla_mask

# TODO: remove `imap` hardcode
"find overlap of IPN/raphe proj field in medulla."
function get_percent_overlap(threshold)
    mm = convert(Array{Bool}, medulla_mask)
    df_f = imap .* mm
    raphe_hits = imap[:,:,:,2] .> threshold
    ipn_hits = imap[:,:,:,1] .> threshold
    sum(ipn_hits .& raphe_hits) / sum(ipn_hits .| raphe_hits)
end

"calc false discovery rate based on threshold compared to control"
function get_fdr(threshold)
    mm = convert(Array{Bool}, medulla_mask)
    df = imap .* mm
    # df = imap
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
z = 5
img = RGB.(zeros(tseriesH, tseriesW))
# channelview(img)[1,:,:,:] .= imadjustintensity(imap[:,:,:,2])
# channelview(img)[2,:,:,:] .= imadjustintensity(imap[:,:,:,1])
im1 = prettify_img(imap[:,:,z,1])
im2 = prettify_img(imap[:,:,z,2])


channelview(img)[1,:,:,:] .= im1 .> reasonable_thresh
channelview(img)[2,:,:,:] .= im2 .> reasonable_thresh
overlap_yrange = 140:380
overlap_xrange = xrange[1]:512-60
img = img[overlap_yrange,overlap_xrange]

img
figW = Float64(uconvert(u"inch", 6u"cm") / u"inch")
figH = figW * length(overlap_yrange) / length(overlap_xrange)
fig, ax = plt.subplots(figsize=(figW,figH))

ax.imshow(permutedims(channelview(img),(2,3,1)),
    interpolation="none")


ax.set_axis_off()
plt.tight_layout()
scalebar3 = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", color="white", location="upper right", font_properties=Dict("size" => 7))
ax.add_artist(scalebar3)
path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_overlap")
@show path*".png"
ax.text(5,length(overlap_yrange)[end]-15,"z = $(Int(etl_vals[z]))μm",fontsize=7, color="white")

fig.savefig(path*".svg", dpi=300)
fig.savefig(path*".png", dpi=300)
fig.savefig(path*".pdf", dpi=300)

fig
# save(path*".png", img)

## save region masks using CMTK
# @pun (fish_dir, region_mask_path, zseries, zbrain_masks) = recording2;
# cmtk_transform_path = joinpath(fish_dir, "warp_out.xform")
# isdir(cmtk_transform_path)
# cmtk_region_masks = L.save_region_masks(region_mask_path, zseries, zbrain_masks, cmtk_transform_path, :cmtk);
