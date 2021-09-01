ENV["DISPLAY"] = "localhost:11.0"
##
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays, AxisArrays
import Lensman: @pun, @assign
using Thunks
import Base.Threads: @spawn, @sync, @threads

L = Lensman

## Load per-region masks from atlas
resources = Resources(zbrain_dir="/data/b115/atlas/zbrain");
@pun (zbrain_dir, zbrain_masks) = resources;

## Need to get `imaging2zseries_plane` (only for visualization)
options = Dict(
    :tseries_dset => nothing, :rel_analysis_dir => "plots",
    :tseries_root_dirs => ["/data/b115"], :slm_root_dirs => ["/data/b115/SLM_files"]
)
recording = Recordings[
        # "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
        # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"

        # "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"  # changing to this one
    ](;options...,
    resources...);
@pun (imaging2zseries_plane) = recording;

## Define directories

# zseries_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/zseries.nrrd"
# transform_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/warp_out.xform"
# region_mask_output_dir = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/h5_output/penis.h5"

# zseries_path = "/scratch/allan/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_fixed.nrrd"
# transform_path = "/scratch/allan/b115/2021-06-08_rsChRmine_h2b6s/fish2/warp_out.xform"
# region_mask_output_dir = "/scratch/allan/b115/2021-06-08_rsChRmine_h2b6s/fish2/h5_output/warped_region_masks.h5"

# zseries_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/20210801T125422241Z_fixed.nrrd"
# transform_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/warp_out.xform"
# region_mask_output_dir = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/h5_output/warped_region_masks.h5"

fish_dir = "/scratch/allan/b115/2021-06-01_rsChRmine_h2b6s/fish3"
zseries_path = joinpath(fish_dir, "zseries.nrrd")
transform_path = joinpath(fish_dir, "warp_out.xform")
region_mask_output_dir = joinpath(fish_dir, "h5_output/warped_region_masks.h5")

# fish_dir = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish2"
# zseries_path = joinpath(fish_dir, "20210801T125422241Z_fixed.nrrd")
# transform_path = joinpath(fish_dir, "warp_out.xform")
# region_mask_output_dir = joinpath(fish_dir, "h5_output/warped_region_masks.h5")

## Test with optic tectum

i = 296  # 106  # index of "Mesencephalon - Tectum Stratum Periventriculare"
name = replace(zbrain_masks["MaskDatabaseNames"][i], "/" => "_")
println("MASK NAME: $name")
mask = read_mask(zbrain_masks, i; rostral=:right, dorsal=:up)
mask = AxisArray(Float32.(mask), AxisArrays.axes(mask))
mask = L.apply_cmtk_transform(zseries_path, mask, transform_path)
mask = mask.data .> 0  # This thing visualizes correctly, so should be good

## Write single mask manually
# h5 = h5open(region_mask_output_dir, "r+")
# mask = sparse(reshape(mask, length(mask), 1))
# using SparseArrays
# H5SparseMatrixCSC(h5, name, mask)

## Visualize
zseries = load(zseries_path)

red = mask.data[:, :, imaging2zseries_plane] .> 0
println(size(red))

v = quantile(zseries[:, :, imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:, :, imaging2zseries_plane], (0, v))
println(size(green))

H, W, nplanes = size(red)
im = zeros(3, H, W, nplanes)

im[1, :, :, :] = red[:, :, :]
im[2, :, :, :] = green[:, :, :]
im = colorview(RGB, im)

##
save("tectum.png", im[:, :, 5])

## Transform all masks

h5 = L.save_region_masks_cmtk(region_mask_output_dir, zseries_path, zbrain_masks, transform_path)

# Closing might fix the issue
close(h5)
