using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed
import PyPlot
import Lensman: @pun, @assign
plt = PyPlot
matplotlib = plt.matplotlib
np = pyimport("numpy")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
# init_workers(18)
##
resources = Resources();
## takes 20 seconds...? should profile macro code...
r = Recordings[
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    ](;resources...);
##
@pun region_mask_path = r
isfile(region_mask_path) ? rm(region_mask_path) : nothing
region_dir = region_mask_path[1:end-3]
isdir(region_dir) ? rm(region_dir, recursive=true) : nothing
# @pun region_masks = r
@pun region_mask_h5 = r
# println("region_masks: $region_masks")
println("success!")
exit()
##
@pun (region_mask_path, zseries, zbrain_masks,
    mm_transform_affine, mm_transform_SyN, zbrain_units) = r;
ants_warps = [mm_transform_affine, mm_transform_SyN];
##
zbrain_masks = masks["MaskDatabaseNames"][1,:]