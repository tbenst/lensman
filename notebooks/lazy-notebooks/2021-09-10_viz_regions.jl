ENV["DISPLAY"] = "localhost:11.0"
##
using AlgebraOfGraphics, CairoMakie
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays, Unitful
using ColorTypes, ColorBlendModes
import Lensman: @pun, @assign
using ImageView
using AxisArrays
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
aog = AlgebraOfGraphics
Data = aog.data

using PyPlot
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

set_aog_theme!() # src
init_workers(16)
##
resources = Resources();
@pun (zbrain_dir, zbrain_masks) = resources;

##
# first fish used in initial submission
# recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
recording_name = "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
tyh5_path = joinpath("/data/dlab/b115", recording_name * "_kalman.h5")
options = Dict(
)
recording = Recordings[recording_name](;options..., resources...);

##
@pun fish_dir = recording

@pun etl_vals = recording;

@pun (imaging2zseries_plane, window_len,
    zbrain_mask_names, fish_dir, zseries, zseries_xml, region_masks_h5
) = recording;
zseries_units = read_microns_per_pixel(zseries_xml);
## add dorsal / ventral 
@pun (region_masks_h5, imaging2zseries_plane, zbrain_mask_names,
    tseries, window_len,
    stim_start_idx, stim_end_idx,
    trial_order, tseriesH, tseriesW, tseriesZ, plot_dir

) = recording;

habenula_name, habenula = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Habenula")
lhemisphere_name, lhemisphere = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Left hemisphere")
rhemisphere_name, rhemisphere = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Right hemisphere")

hab = read_registered_mask(region_masks_h5, habenula_name);
left_hemi = read_registered_mask(region_masks_h5, "Left hemisphere");
right_hemi = read_registered_mask(region_masks_h5, "Right hemisphere");

    
prof = sum(hab,dims=[1,2])[1,1,:]
cumprof = cumsum(prof)
halfway = searchsortedfirst(cumprof, cumprof[end] / 2)


# ymax = maximum(prof)
# df2 = (z=[halfway, halfway], area=[0, ymax])
# fig = Figure()
# ax = Makie.Axis(fig[1, 1], title="dorsal/ventral hab cutoff")
# p1 = (aog.data((z=1:length(prof), area=prof)) + aog.data(df2) ) * mapping(:z, :area)* visual(Lines)
# aog.draw!(fig[1,1], p1)
# fig
dorsal_hab = copy(hab)
ventral_hab = copy(hab)
dorsal_hab[:,:,halfway+1:end] .= false
ventral_hab[:,:,1:halfway] .= false
# im = RGB.(dorsal_hab)
# channelview(im)[[2,3],:,:,:] .= 0
# channelview(im)[2,:,:,:] = ventral_hab
# imshow(im)
masks = [dorsal_hab .& left_hemi, dorsal_hab .& right_hemi,
    ventral_hab .& left_hemi, ventral_hab .& right_hemi]
masks = map(m->imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);
region_names = ["left Diencephalon - Dorsal Habenula", "right Diencephalon - Dorsal Habenula",
    "left Diencephalon - Ventral Habenula", "right Diencephalon - Ventral Habenula"]

##
# choose proper z-plane (look at titration notebook)
z = 3
zseriesZplane = imaging2zseries_plane[z]
# list of regions to show (check out lensman code)
##
# order is critical
REGION_LIST = [
    "Glyt2 Stripe 3"
    # "Cerebellum"
    "Raphe - Superior"
    "Rhombencephalon - Gad1b Cluster 16"
    "Rhombencephalon - Gad1b Cluster 2"
    # "Rhombencephalon - Gad1b Cluster 1"
    # "Glyt2 Cluster 12"
    "Tectum Stratum Periventriculare"
    # "Rhombomere 1"
    # "Subpallial Gad1b cluster"
]
region_masks = [L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane, region) for region in REGION_LIST];
insert!(region_masks, 1, ("Dorsal Habenula", dorsal_hab[:,:,imaging2zseries_plane]))
insert!(region_masks, 2, ("Ventral Habenula", ventral_hab[:,:,imaging2zseries_plane]))
# @assert length(region_masks) == 11

# colormap choose seperate colors
# imshow zseries (matplotlib), add scalebar
# add color to zseries image, by taking maximum value..?

##
v = quantile(zseries[:,:,imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:,:,imaging2zseries_plane], (0, v))
thresh = otsu_threshold(green[:]) * 0.5
for c in findall(green .<= thresh)
    green[c] = 0
end
# green = adjust_histogram(green,GammaCorrection(0.5))
im2 = RGB.(green)
im2
##

seeds = vcat([RGB(0,x,0) for x in 0:0.11:1]..., RGB(1,1,1), RGB(0,0,0))
region_colors = distinguishable_colors(length(region_masks), seeds, dropseed=true)

region_masks

im2 = RGB.(green)
channelview(im2)[[1,3],:,:,:] .= 0
for (r, reg) in enumerate(region_masks)
    # region = reg .> 0

    color_mask = RGB{N0f16}.(zeros(size(im2)...))
    color_mask[reg[2]] .= region_colors[r]

    im2 = mapc.((a,b)->maximum([a,b]), im2, color_mask)
    im2 = mapc.((a,b)->b==0 ? a : b, im2, color_mask)

    # im2[reg[2]] .= region_colors[r]
    im2 = blend.(im2, color_mask, mode=BlendScreen)
    # im2 = blend.(im2, color_mask, mode=BlendLighten)
    # im2 = blend.(im2, color_mask, mode=BlendColor)
end
im2[:,:,z]
# imshow(im2)

#
zs = [4, 6, 8, 10]
# no idea why need to scale up
figB = Float64(uconvert(u"inch", 55u"mm")) / 1u"inch" * 1.3
fig, axs = plt.subplots(2,2,figsize=(figB,0.95*figB))
for (z,ax) in zip(zs,np.transpose(axs))
    py_im = permutedims(channelview(im2[:,:,z]),[2,3,1])
    ax.imshow(py_im, interpolation="none")
    ax.set_axis_off()
end
plt.tight_layout()
fig.subplots_adjust(wspace=0.05, hspace=0.05)

um_per_px1 = zseries_units[2] / 1u"μm"
scalebar1 = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="lower right", color = "white",
    font_properties=Dict("size" => 7))
axs[1].add_artist(scalebar1)
plotpath = joinpath(plot_dir,"region_annotations")
@show plotpath*".svg"
fig.savefig(plotpath*".svg",dpi=300)
fig.savefig(plotpath*".pdf",dpi=300)
fig.savefig(plotpath*".png",dpi=300)
fig
##
region_colors
##
@show etl_vals[zs]

## assess drift

img = RGB.(adjust_histogram(
    imadjustintensity(zseries[:,:,imaging2zseries_plane[10]]),
    GammaCorrection(0.2)))
mean_tseries = mean(tseries[:,:,10,1:10:1000],dims=3)[:,:,1];
channelview(img)[3,:,:] .= 0
channelview(img)[1,:,:] .= adjust_histogram(
    Gray.(imadjustintensity(imresize(mean_tseries,1024,1024),
    )), GammaCorrection(0.2))
img