ENV["DISPLAY"] = "localhost:11.0"
##
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays
import Lensman: @pun, @assign
using ImageMagick, Images, ImageView
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman

##
resources = Resources(zbrain_dir="/data/b115/atlas/zbrain");
@pun (zbrain_dir, zbrain_masks) = resources;

##
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

##
# file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/zseries.nrrd"
# zseries_path = "/data/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_fixed.nrrd"  # changing to this one
# zseries_path = "/data/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_fixed.nrrd"  # changing to this one

zseries_path = "/data/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/20210801T125422241Z_fixed.nrrd"
zseries = load(zseries_path)

# Load atlas
# file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/zbrain_h2b.nrrd"
# zseries_h2b = load(file_path)

##

# Load warped atlas (baseline is ants registration)
# using NIfTI
# file_path = "/data/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/20210726T095924956Z_Warped.nii.gz"
# zseries_ants = niread(file_path)

## Load warped atlas
# file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/warped_boy.nrrd"
# file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/warped_full_new.nrrd"
# file_path = "/scratch/allan/b115/2021-06-08_rsChRmine_h2b6s/fish2/warped.nrrd"  # changing to this one

file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/warped.nrrd"
zseries_cmtk = load(file_path)

println(minimum(zseries_cmtk))
println(maximum(zseries_cmtk))

##
# Visualize warped (now that at least not all 0's)
# v = quantile(zseries_cmtk[:, :, imaging2zseries_plane][:], 0.99)
# red = imadjustintensity(zseries_cmtk[:, :, imaging2zseries_plane], (0, v))

# println(size(zseries_cmtk))
# im2 = zseries_cmtk[:, :, imaging2zseries_plane]
# save("cmtk_affine.png", im2[:, :, 5])

##
zseries_cmtk_data = convert(Array{Float32,3}, zseries_cmtk.data)
v = quantile(zseries_cmtk_data[:, :, imaging2zseries_plane][:], 0.99)
println(v)
red = imadjustintensity(zseries_cmtk_data[:, :, imaging2zseries_plane], (0, v))  # can use its own quantile
println(size(red))

v = quantile(zseries[:, :, imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:, :, imaging2zseries_plane], (0, v))
println(size(green))

H, W, nplanes = size(red)
im = zeros(3, H, W, nplanes)

im[1, :, :, :] = red[:, :, :]
im[2, :, :, :] = green[:, :, :]
im = colorview(RGB, im)
save("warped_full_new_fish2.png", im[:, :, 5])

##
v = quantile(zseries_ants[:, :, imaging2zseries_plane][:], 0.99)
println(v)
red = imadjustintensity(zseries_ants[:, :, imaging2zseries_plane], (0, v))  # can use its own quantile
println(size(red))

v = quantile(zseries[:, :, imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:, :, imaging2zseries_plane], (0, v))
println(size(green))

H, W, nplanes = size(red)
im = zeros(3, H, W, nplanes)

im[1, :, :, :] = red[:, :, :]
im[2, :, :, :] = green[:, :, :]
im = colorview(RGB, im)
save("penis2.png", im[:, :, 5])

##

# v = quantile(zseries_ants[:, :, imaging2zseries_plane][:], 0.99)
# green = imadjustintensity(zseries_ants[:, :, imaging2zseries_plane], (0, v))
# println(size(green))

H, W, nplanes = size(red)
im = zeros(3, H, W, nplanes)

im[1, :, :, :] = red[:, :, :]
im[2, :, :, :] = green[:, :, :]
im = colorview(RGB, im)
save("penis.png", im[:, :, 5])

##

v = quantile(zseries_ants[:, :, imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries_ants[:, :, imaging2zseries_plane], (0, v))

# im2 = RGB.(green)
# channelview(im2)[1,:,:,:] .= red
# save("penis.png", im2[:, :, 5])

##
imshow(im2)




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
if mean(target_groups[1], dims=1)[2] > 256
    stim_names = ["right habenula", "left habenula", "control"]
else
    stim_names = ["left habenula", "right habenula", "control"]
end
region_names = ["habenula", "optic tectum", "rhombomere 1", "raphe", "cerebellum"]
masks = [habenula, tectum, rhom1, raphe, cerebellum];
masks = vcat(map(m -> [m .& rhemisphere, m .& lhemisphere], masks)...)
region_names = vcat(map(n -> ["right $n", "left $n"], region_names)...)
# region_names = vcat(map(n->[n, n], region_names)...)
masks = map(m -> imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);
# imshow(masks[1]) # right hab works
##

# @everywhere include("/home/tyler/code/lensman/src/tidying.jl")
# @everywhere using DataFrames, Statistics, HDF5
regions_df = L.per_trial_regions_df(tseries, window_len,
    stim_start_idx, stim_end_idx,
    # stimStartIdx, stimEndIdx,
    trial_order, masks, region_names, stim_names)

##
fig = Figure(resolution=(1200, 500))
p1 = Data(regions_df) * visual(BoxPlot, show_notch=true) *
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