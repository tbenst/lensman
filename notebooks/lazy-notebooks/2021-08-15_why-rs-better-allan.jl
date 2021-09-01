ENV["DISPLAY"] = "localhost:13.0"
##
using AlgebraOfGraphics, CairoMakie
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays
import Lensman: @pun, @assign
using ImageView
using AxisArrays
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
aog = AlgebraOfGraphics
Data = aog.data
set_aog_theme!() # src
init_workers(16)

##
# The good one
recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"

# recording_name = "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"

##
resources = Resources(zbrain_dir="/data/b115/atlas/zbrain");
@pun (zbrain_dir, zbrain_masks) = resources;

##
# NOTE: need to make sure to have copied SLM files from /mnt/deissero to /data/b115/SLM_files/14-Jul-2021
# at least user=allan can't read them from /mnt/deissero
tyh5_path = "/data/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5"
options = Dict(
    # :rel_analysis_dir => "plots", :h5_read_strategy => :lazy_ty5, :lazy_tiff => false,  # just in case of trouble later: added `:lazy_tiff => false` here
    :rel_analysis_dir => "plots-kalman", :h5_read_strategy => :lazy_hwzt, :tyh5_path => tyh5_path, :tseries_dset => "kalman",
    :tseries_root_dirs => ["/data/b115"], :slm_root_dirs => ["/data/b115/SLM_files"]
)
recording = Recordings[recording_name](;options..., resources...);

##
@pun fish_dir = recording

@pun (tseries, tseriesH) = recording;

      @pun (imaging2zseries_plane, window_len, 
      tseries, tseriesH, tseriesW, tseriesZ, trial_order, 
      zbrain_mask_names, target_groups, fish_dir,
      stim_start_idx, stim_end_idx) = recording;

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

## 
# Remove masks returning as `nothing` (assuming correct flushing when creation region_masks.h5, this should
# just be masks saved as 0 because no pixels overlapping with imaged - should only be a few)
valid_mask_indices = (masks .!== nothing)
println("Dropping $(length(masks) - sum(valid_mask_indices)) valid masks out of $(length(masks)), because = nothing")
region_names = region_names[valid_mask_indices];
masks = masks[valid_mask_indices];

##
# Remove masks with no voxels imaged
valid_mask_indices = map(m -> sum(m) > 0, masks)
println("Dropping a further $(length(masks) - sum(valid_mask_indices)) valid masks out of $(length(masks)), because 0 in imaged plane")
region_names = region_names[valid_mask_indices];
masks = masks[valid_mask_indices];

##
# Read hemisphere_masks
lhemisphere_name, lhemisphere = L.read_first_mask(
    region_masks_h5, zbrain_mask_names, imaging2zseries_plane, "Left hemisphere");
rhemisphere_name, rhemisphere = L.read_first_mask(
    region_masks_h5, zbrain_mask_names, imaging2zseries_plane, "Right hemisphere");

##
# Construct "left" and "right" masks
masks = vcat(map(m -> [m .& rhemisphere, m .& lhemisphere], masks)...);
region_names = vcat(map(n -> ["right $n", "left $n"], region_names)...);

# Resize masks to tseries shape
masks = map(m -> imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);

##
# FIXME: this right now is setup to run with tyh5, minor change for tiff's
regions_df = L.per_trial_regions_df(
    tseries,
    window_len, stim_start_idx, stim_end_idx,
    trial_order, masks, region_names)  # , stim_names=nothing)
# regions_df = L.per_trial_regions_df(
#     tseries,
#     window_len, stim_start_idx, stim_end_idx,
#     trial_order, masks, region_names)  # , stim_names=nothing)

##

# Path to write out to (in the case that tyh5_path is defined, otherwise reading from tiff's)
out_dir = joinpath(mask_dir, basename(recording[:tyh5_path]), basename(recording[:tseries_dset]))
# out_dir = joinpath(mask_dir, basename(recording_name), "raw")
println("Will create output_dir='$out_dir'")

##
mkpath(out_dir)

##
using Arrow
open(joinpath(out_dir, "regions_df.arrow"), "w") do io
    Arrow.write(io, regions_df)
end


#### RUNNING UNTIL HERE IS GOOD FOR CREATING DATAFRAME

#########################################################################################################

##
tectum_name, tectum = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Tectum")
cerebellum_names, cerebellum = L.read_matching_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Cerebe")
medulla_names, medulla = L.read_matching_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Medulla")
raphe_names, raphe = L.read_matching_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Raphe")


habenula_name, habenula = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Habenula")
rhom1_name, rhom1 = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombomere 1")
lhemisphere_name, lhemisphere = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Left hemisphere")
rhemisphere_name, rhemisphere = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Right hemisphere")
tectum_name, medulla_names, cerebellum_names, raphe_names, rhom1_name, habenula_name

##
if mean(target_groups[1], dims=1)[2] > 256
    stim_names = ["right habenula", "left habenula", "control"]
else
    stim_names = ["left habenula", "right habenula", "control"]
end
region_names = ["habenula", "optic tectum", "rhombomere 1", "raphe", "cerebellum"]
masks = [habenula, tectum, rhom1, raphe, cerebellum];
masks = vcat(map(m -> [m .& rhemisphere, m .& lhemisphere], masks)...)
region_names = vcat(map(n -> ["right $n", "left $n"], region_names)...)
masks = map(m -> imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);

##

# regions_df = L.per_trial_regions_df(tseries, window_len,
#     stim_start_idx, stim_end_idx,
#     trial_order, masks, region_names)  # , stim_names=nothing)            # CAN CHANGE BACK FROM THIS

regions_df = L.per_trial_regions_df(
    recording[:tyh5_path], recording[:tseries_dset],
    window_len, stim_start_idx, stim_end_idx,
    trial_order, masks, region_names)  # , stim_names=nothing)

## Save DataFrame (without filtering any possible NaN's for now)

# Path to write out to (in the case that tyh5_path is defined, otherwise reading from tiff's)
out_dir = joinpath(mask_dir, basename(recording[:tyh5_path]), basename(recording[:tseries_dset]))
println("Will create output_dir='$out_dir'")

##
mkdir(out_dir)

##
open(joinpath(out_dir, "regions_df_tiny.arrow"), "w") do io
    Arrow.write(io, regions_df)
end

##
regions_df = filter(:("Δf/f") => !isnan, regions_df)

##
fig = Figure(resolution=(1200, 500))
# p1 = Data(regions_df) * visual(BoxPlot, show_notch=true) *
#     mapping(:region, "Δf/f", color=:hemisphere, dodge=:hemisphere, col=:stim)
p1 = Data(filter(row -> row.stim == "1" || row.stim == "16", regions_df)) * visual(BoxPlot, show_notch=true) *
    mapping(:region, "Δf/f", color=:hemisphere, dodge=:hemisphere, col=:stim)
grid = draw!(fig, p1, axis=(xticklabelrotation = pi / 6,))
legend!(fig[1,end + 1], grid)
ppath = joinpath(plot_dir, "region_df_f_by-trial")
save(ppath * ".png", fig)
save(ppath * ".svg", fig)
fig
##
regions_df[:,:fish] .= recording_folder * fish_name
regions_df[:,:experiment] .= exp_name
regions_df[:,:uri] .= recording[:uri]
open(df_dir * "regions_df.arrow", "w") do io
    Arrow.write(io, regions_df)
end
##
v = quantile(zseries[:,:,imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:,:,imaging2zseries_plane], (0, v))
thresh = otsu_threshold(green[:]) * 0.5
for c in findall(green .<= thresh)
    green[c] = 0
end
# green = adjust_histogram(green,GammaCorrection(0.5))
im2 = RGB.(green)
channelview(im2)[1,:,:,:] .= habenula .* 0.8
channelview(im2)[1,:,:,:] .+= rhom1 .* 0.8
channelview(im2)[3,:,:,:] .= tectum .* 0.8
# imshow(im2)
im2[:,:,6]
save(joinpath(plot_dir, "tectum_rhomb1.png"), im2[:,:,6])
##
im3 = RGB.(green)
channelview(im3)[1,:,:,:] .= lhemisphere .* 0.5
channelview(im3)[3,:,:,:] .= rhemisphere .* 0.5
save(joinpath(plot_dir, "hemishpere.png"), im3[:,:,6])
im3[:,:,6]