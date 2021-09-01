ENV["DISPLAY"] = "localhost:11.0"
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

##
recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"

##
resources = Resources(zbrain_dir="/data/b115/atlas/zbrain");
@pun (zbrain_dir, zbrain_masks) = resources;

##
# NOTE: need to make sure to have copied SLM files from /mnt/deissero to /data/b115/SLM_files/14-Jul-2021
# at least user=allan can't read them from /mnt/deissero
options = Dict(
    # FIXME: no bueno
    :tseries_dset => nothing, :rel_analysis_dir => "plots", :lazy_tiff => false,  # use raw data
    :tseries_root_dirs => ["/data/b115"], :slm_root_dirs => ["/data/b115/SLM_files"]
)
recording = Recordings[recording_name](;options..., resources...);

##
@pun (imaging2zseries_plane, window_len, 
      tseries, tseriesH, tseriesW, tseriesZ, trial_order, 
      zbrain_mask_names, target_groups,
      stim_start_idx, stim_end_idx) = recording;

## Open h5 for warped region masks (saved to standardized path)
scratch_dir = "/scratch/allan/b115"
mask_dir = joinpath(scratch_dir, dirname(recording_name), "h5_output")
println("Will read masks from: '$mask_dir'")
region_masks_h5_path = joinpath(mask_dir, "warped_region_masks.h5")
region_masks_h5 = h5open(region_masks_h5_path, "r", swmr=true)
plot_dir = joinpath(mask_dir, "plots")
println("Will save plots to '$plot_dir'")

##

# Read each mask and iteratively call `per_trial_regions_df`, and vcat the results into one dataframe. 
# Then can do something like:
# `median_region_df = combine(groupby(regions_df, [:region, :stim, :hemisphere]), "Δf/f" => median)`

# calc_diff_between_lrhab_and_control_stim(df) = maximum(df[df[:,:stimNum].==2] - df[df[:,:stimNum].==3], df[df[:,:stimNum].==1] - df[df[:,:stimNum].==3])

# stim_effectiveness[argsort(stim_effectiveness[:,"Δf/f_median" ],:]

function read_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, mask_name)
    _, mask = L.read_first_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, mask_name)
    mask
end

# TODO: Remove RH and LH (can simply go up to end - 2)
region_names = zbrain_masks["MaskDatabaseNames"][1:end - 2]
masks = map(n -> read_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, n),
            region_names)

##
valid_mask_indices = (masks .!== nothing)
println("There are $(sum(valid_mask_indices)) valid masks")
region_names = region_names[valid_mask_indices]
masks = masks[valid_mask_indices]
masks = map(m -> imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);

# Right now 23 are nothing (this could just be a matter of re-running with close(h5) to make sure
# flushed, but more likely some masks are just empty)
# NOTE: for now not doing "left" and "right", but instead just getting coarse sense of most important
# regions

##
# regions_df = L.per_trial_regions_df(tseries, window_len,
#     stim_start_idx, stim_end_idx,
#     trial_order, masks, region_names)  # , stim_names=nothing)
regions_df = L.per_trial_regions_df(
    recording[:tyh5_path], recording[:tseries_dset],
    window_len, stim_start_idx, stim_end_idx,
    trial_order, masks, region_names)  # , stim_names=nothing)

## Write to Arrow
using Arrow
open(joinpath(plot_dir, "regions_df_1_16_smoothed.arrow"), "w") do io  # make name depend on whether data was smoothed or not
    Arrow.write(io, regions_df)
end

##

regions_df = filter(row -> row.stim == "1" || row.stim == "16", regions_df)

## Drop NaN's for cleanliness
regions_df = filter(:("Δf/f") => !isnan, regions_df)


##
median_region_df = combine(groupby(regions_df, [:region, :stim]), "Δf/f" => median)

##
df_baseline = filter(row -> row.stim == "1", median_region_df)
df = filter(row -> row.stim == "16", median_region_df)
df[:, :effectiveness] = df[:,"Δf/f_median"] .- df_baseline[:,"Δf/f_median"]

##
sort!(df, [:effectiveness])

##
function effectiveness(df)
    median(df[df[:, :stim] .== 16] - df[df[:, :stim] .== 1]) 
end

# stim_effectiveness = map(effectiveness, groupby(median_region_df, :region))
stim_df = combine(groupby(median_region_df, :region), "Δf/f_median" .=> [maximum minimum])

# stim_effectiveness[argsort(stim_effectiveness[:,"Δf/f_median" ],:]
##
stim_df[:horse] = stim_df[:, "Δf/f_median_maximum"] - stim_df[:, "Δf/f_median_minimum"]

# Getting there