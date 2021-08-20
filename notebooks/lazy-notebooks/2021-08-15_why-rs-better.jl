ENV["DISPLAY"] = "localhost:11"
##
using AlgebraOfGraphics, CairoMakie
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays
import Lensman: @pun, @assign
using ImageView
using AxisArrays
# import PyPlot
# plt = PyPlot
# matplotlib = plt.matplotlib
# np = pyimport("numpy")
# import Plots
# import Plots: heatmap
# using Gadfly
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
aog = AlgebraOfGraphics
Data = aog.data
set_aog_theme!() #src
init_workers(36)


##
resources = Resources();
##
@pun (zbrain_dir, zbrain_masks) = resources
options = Dict(
    :tseries_dset=>nothing, :rel_analysis_dir=>"plots" # use raw data
)
recording = Recordings[
        "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
    ](;options...,
    resources...);
recording.nodes[:zbrain_registered].restore
recording.nodes[:zbrain_registered].checkpoint
# recording.nodes[:zbrain_registered].wrapped_thunk
##
@pun region_masks_h5 = recording
# @pun (nCells, cell_centers, cells_mask, iscell) = r
##
@pun (region_masks_h5, imaging2zseries_plane,
    zbrain_masks, zbrain_mask_names, plot_dir, 
    tseries, trial_average, window_len, trial_order,
    stim_start_idx, stim_end_idx, target_groups,
    tyh5_path, tseries_dset, zbrain_transforms,
    tseriesH, tseriesW, tseriesZ, zseries, rostral, dorsal,
    recording_folder, fish_name,
    exp_name
) = recording;
##
imap = influence_map(trial_average, window_len);
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

# stim_names = ["left habenula", "right habenula", "control"]
if mean(target_groups[1],dims=1)[2] > 256
    stim_names = ["right habenula", "left habenula", "control"]
else
    stim_names = ["left habenula", "right habenula", "control"]
end
region_names = ["habenula", "optic tectum", "rhombomere 1", "raphe", "cerebellum"]
masks = [habenula, tectum, rhom1, raphe, cerebellum];
masks = vcat(map(m->[m .& rhemisphere, m .& lhemisphere], masks)...)
region_names = vcat(map(n->["right $n", "left $n"], region_names)...)
# region_names = vcat(map(n->[n, n], region_names)...)
masks = map(m->imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);
# imshow(masks[1]) # right hab works
##

# @everywhere include("/home/tyler/code/lensman/src/tidying.jl")
# @everywhere using DataFrames, Statistics, HDF5
regions_df = L.per_trial_regions_df(tseries, window_len,
    stim_start_idx, stim_end_idx,
    # stimStartIdx, stimEndIdx,
    trial_order, masks, region_names, stim_names)

##
fig = Figure(resolution = (1200, 500))
p1 = Data(regions_df) * visual(BoxPlot, show_notch=true) *
    mapping(:region, "Δf/f", color=:hemisphere, dodge=:hemisphere, col=:stim)
grid = draw!(fig,p1, axis=(xticklabelrotation = pi/6,))
legend!(fig[1,end+1], grid)
ppath = joinpath(plot_dir, "region_df_f_by-trial")
save(ppath*".png", fig)
save(ppath*".svg", fig)
fig
##
regions_df[:,:fish] .= recording_folder * fish_name
regions_df[:,:experiment] .= exp_name
regions_df[:,:uri] .= recording[:uri]
open(df_dir*"regions_df.arrow", "w") do io
    Arrow.write(io, regions_df)
end
##
v = quantile(zseries[:,:,imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:,:,imaging2zseries_plane], (0,v))
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