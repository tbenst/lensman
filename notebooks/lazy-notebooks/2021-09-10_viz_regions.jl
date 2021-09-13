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
resources = Resources(zbrain_dir="/data/b115/atlas/zbrain");
@pun (zbrain_dir, zbrain_masks) = resources;

##
# DONE
# recording_name = "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
# recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-118trial-122"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"

# DONE - but errors (masks are same across planes; only 23 nonzero)
# recording_name = "2021-06-15_rsChRmine_h2b6s_6dpf/fish1/TSeries-titration-192trial-050"

# can't do this one as forgot to acquire ZSeries
# recording_name = "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061"
##
# NOTE: need to make sure to have copied SLM files from /mnt/deissero to /data/b115/SLM_files/14-Jul-2021
# at least user=allan can't read them from /mnt/deissero
# tyh5_path = "/data/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5"
tyh5_path = joinpath("/data/b115", recording_name * "_kalman.h5")
options = Dict(
    # :rel_analysis_dir => "plots", :h5_read_strategy => :lazy_ty5, :lazy_tiff => false,  # just in case of trouble later: added `:lazy_tiff => false` here
    :rel_analysis_dir => "plots-kalman", :h5_read_strategy => :lazy_hwzt, :tyh5_path => tyh5_path, :tseries_dset => "kalman",
    :tseries_root_dirs => ["/data/b115"],
    :slm_root_dirs => ["/data/b115/SLM_files"]
)
recording = Recordings[recording_name](;options..., resources...);

##
@pun fish_dir = recording

@pun (tseries, tseriesH, etl_vals) = recording;

@pun (imaging2zseries_plane, window_len, tseriesH, tseriesW, tseriesZ,
    zbrain_mask_names, fish_dir, zseries, zseries_xml
) = recording;
zseries_units = read_microns_per_pixel(zseries_xml);
## Open h5 for warped region masks (saved to standardized path)

scratch_dir = "/scratch/allan/b115"
mask_dir = joinpath(scratch_dir, dirname(recording_name), "h5_output")
println("Will read masks from: '$mask_dir'")
region_masks_h5_path = joinpath(mask_dir, "warped_region_masks.h5")
region_masks_h5 = h5open(region_masks_h5_path, "r", swmr=true)
plot_dir = joinpath(mask_dir, "plots")
println("Will save plots to '$plot_dir'")
#########################################################################################################
##
# This function is basically nothing (I was planning on doing something more complicated initially)
function read_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, mask_name)
    _, mask = L.read_first_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, mask_name)
    mask
end

# Remove RH and LH (can simply go up to end - 2)
region_names = zbrain_masks["MaskDatabaseNames"][1:end - 2];
masks = map(n -> read_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, n),
            region_names);
## TODO: 
# choose proper z-plane (look at titration notebook)
z = 3
zseriesZplane = imaging2zseries_plane[z]
# list of regions to show (check out lensman code)
##
# order is critical
map(r->r[1], region_masks)
REGION_LIST = [
    "Cerebellum"
    "Rhombencephalon - Gad1b Cluster 1"
    "Rhombencephalon - Gad1b Cluster 16"
    "Rhombencephalon - Gad1b Cluster 2"
    "Glyt2 Cluster 12"
    "Glyt2 Stripe 3"
    "Habenula"
    "Tectum Stratum Periventriculare"
    "Raphe - Superior"
    "Rhombomere 1"
    "Subpallial Gad1b cluster"
]
region_masks = [L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane, region) for region in REGION_LIST];

# colormap choose seperate colors
# imshow zseries (matplotlib), add scalebar
# add color to zseries image, by taking maximum value..?

##


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
colors = distinguishable_colors(length(region_masks), seeds, dropseed=true)

region_masks

im2 = RGB.(green)
channelview(im2)[[1,3],:,:,:] .= 0
for (r, reg) in enumerate(region_masks)
    # region = reg .> 0
    color_mask = RGB{N0f16}.(zeros(size(im2)...))
    color_mask[reg[2]] .= colors[r]
    # im2 = mapc.((a,b)->maximum([a,b]), im2, color_mask)
    im2 = blend.(im2, color_mask, mode=BlendScreen, opacity=1.0)
end
im2[:,:,z]
imshow(im2)

##
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

# colors
##
@show etl_vals[zs]