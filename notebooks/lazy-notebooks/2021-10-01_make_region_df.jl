# copied / adopted from 2021-08-15_why-rs-better-allan.jl on Osprey
# we do manual alignment for 2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072
ENV["DISPLAY"] = "localhost:12.0"
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
init_workers(36)
##
resources = Resources();
@pun (zbrain_dir, zbrain_masks) = resources;

##
# DONE
# recording_name = "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072" # wrong date / name, oops
recording_name = "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
# recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
# recording_name = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-118trial-122"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
# recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"

tyh5_path = joinpath("/data/b115", recording_name * "_kalman.h5")
options = Dict()
recording = Recordings[recording_name](;options..., resources...);

@pun tyh5_path = recording
##
@pun fish_dir = recording

@pun (tseries, tseriesH) = recording;

@pun (imaging2zseries_plane, window_len, 
tseries, tseriesH, tseriesW, tseriesZ, trial_order, 
zbrain_mask_names, target_groups, fish_dir,
stim_start_idx, stim_end_idx, nStimuli, zseries) = recording;

@pun (plot_dir, region_masks_h5, cell_masks,cells) = recording;


## manual registration of zseries to tseries
mean_tseries = mean(tseries[:,:,:,1:10:1000], dims=4)[:,:,:,1];
Gray.(imadjustintensity(mean_tseries[:,:,10]))
## assess drift
zplane_offset = Dict{Int,NamedTuple}()
zplane_offset[1] = (y_offset=-15,x_offset=0,z_offset=0)
zplane_offset[2] = (y_offset=-12,x_offset=0,z_offset=0)
zplane_offset[3] = (y_offset=-9,x_offset=0,z_offset=0)
zplane_offset[4] = (y_offset=-6,x_offset=0,z_offset=0)
zplane_offset[5] = (y_offset=5,x_offset=0,z_offset=0)
zplane_offset[6] = (y_offset=10,x_offset=0,z_offset=1)
zplane_offset[7] = (y_offset=15,x_offset=4,z_offset=1)
zplane_offset[8] = (y_offset=20,x_offset=8,z_offset=1)
zplane_offset[9] = (y_offset=26,x_offset=8,z_offset=1)
zplane_offset[10] = (y_offset=30,x_offset=8,z_offset=1)
# cycle through all 10 planes..
z = 10 

shifted_zseries = circshift(zseries, zplane_offset[z])
img = RGB.(adjust_histogram(
    imadjustintensity(shifted_zseries[:,:,imaging2zseries_plane[z]]),
    GammaCorrection(0.2)))

mean_tseries = mean(tseries[:,:,z,1:10:1000],dims=3)[:,:,1];
channelview(img)[3,:,:] .= 0
channelview(img)[1,:,:] .= adjust_histogram(
    Gray.(imadjustintensity(imresize(mean_tseries,1024,1024),
    )), GammaCorrection(0.2))
# channelview(img)[1,:,:] .= 0
# channelview(img)[2,:,:] .= 0
img

##
# Remove RH and LH (can simply go up to end - 2)
region_names = zbrain_masks["MaskDatabaseNames"][1:end - 2];
masks = @showprogress map(n -> L.try_read_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, n),
            region_names);


## add dorsal/ventral habenula
habenula_name, habenula = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Habenula")

hab = read_registered_mask(region_masks_h5, habenula_name);
prof = sum(hab,dims=[1,2])[1,1,:]
cumprof = cumsum(prof)
halfway = searchsortedfirst(cumprof, cumprof[end] / 2)
dorsal_hab = copy(hab)
ventral_hab = copy(hab)
dorsal_hab[:,:,halfway+1:end] .= false
ventral_hab[:,:,1:halfway] .= false
hab_masks = [dorsal_hab, ventral_hab]
hab_masks = map(m->imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, hab_masks);
hab_region_names = ["Diencephalon - Dorsal Habenula", "Diencephalon - Ventral Habenula",]
length(masks)
masks = vcat(masks, hab_masks);
region_names = vcat(region_names, hab_region_names)
## 
# Remove masks returning as `nothing` (assuming correct flushing when creation region_masks.h5, this should
# just be masks saved as 0 because no pixels overlapping with imaged - should only be a few)
valid_mask_indices = (~).(isnothing.(masks))
println("Dropping $(length(masks) - sum(valid_mask_indices)) valid masks out of $(length(masks)), because = nothing")
region_names = region_names[valid_mask_indices];
masks = masks[valid_mask_indices];


# Remove masks with no voxels imaged
valid_mask_indices = map(m -> sum(m) > 0, masks)
println("Dropping a further $(length(masks) - sum(valid_mask_indices)) valid masks out of $(length(masks)), because 0 in imaged plane")
region_names = region_names[valid_mask_indices];
masks = masks[valid_mask_indices];


# Read hemisphere_masks
lhemisphere_name, lhemisphere = L.read_first_mask(
    region_masks_h5, zbrain_mask_names, imaging2zseries_plane, "Left hemisphere");
rhemisphere_name, rhemisphere = L.read_first_mask(
    region_masks_h5, zbrain_mask_names, imaging2zseries_plane, "Right hemisphere");
rhemisphere = imresize(rhemisphere, (tseriesH, tseriesW)) .> 0;
lhemisphere = imresize(lhemisphere, (tseriesH, tseriesW)) .> 0;

## Construct "left" and "right" masks
# Resize masks to tseries shape
masks = map(m -> imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);
masks = vcat(map(m -> [m .& rhemisphere, m .& lhemisphere], masks)...);
region_names = vcat(map(n -> ["right $n", "left $n"], region_names)...);

## make stim masks for exclusion zones...
cell_masks
stim_groups_df = unique(cells, [:cellID, :stimNum])
stim_exclude_masks = zeros(Bool, tseriesH, tseriesW, tseriesZ, nStimuli);

for s in 1:nStimuli
    cellIDs = filter(r->r.stimNum .== s, stim_groups_df).cellID
    for c in cellIDs
        view(stim_exclude_masks,:,:,:,s)[cell_masks[c]] .= true
    end
end;
zproj = maximum(stim_exclude_masks, dims=3)
stim_exclude_masks .= zproj
Gray.(stim_exclude_masks[:,:,5,2])

## viz a mask
idx = findall(occursin.("Raphe", region_names))[end-1]
@show region_names[idx]
# imshow(masks[idx])
mask = masks[idx]
Gray.(mask[:,:,10])
##
mask = L.exclude_from_masks([mask], stim_exclude_masks[:,:,:,2])[1]
Gray.(mask[:,:,10])
# Gray.(masks[idx][:,:,1])

## use manual zseries alignment
# we ignore z offset for now
aligned_masks = deepcopy(masks)
for i in 1:length(aligned_masks)
    for z in 1:size(aligned_masks[i],3)
        (y_offset,x_offset,z_offset) = zplane_offset[z]
        y_offset = Int(floor(y_offset / 2)) # 1024 offset to 512
        view(aligned_masks[i],:,:,z) .=
            circshift(aligned_masks[i][:,:,z], [y_offset,x_offset])
    end
end;
##

per_stim_region_masks = repeat(aligned_masks, 1, nStimuli)
for s in 1:nStimuli
    per_stim_region_masks[:,s] = L.exclude_from_masks(aligned_masks, stim_exclude_masks[:,:,:,s])
end

# sanity check--show 2nd order all-region mask..
Gray.((reduce((a,b)->a .| b, per_stim_region_masks[:,1])))[:,:,1]

rap_idx = findall(occursin.("Raphe", region_names))
rap_names = region_names[rap_idx]
gad_idx = findall(occursin.("Gad1b Cluster 16", region_names))
gad_names = region_names[gad_idx]
@assert length(rap_idx) == 4
hab_idx = findall(occursin.("Habenula", region_names))
hab_names = region_names[hab_idx]
@assert length(hab_idx) == 10
Gray.((reduce((a,b)->a .| b, per_stim_region_masks[:,1])))[:,:,1]
og_raphe = aligned_masks[rap_idx[3]] .| aligned_masks[rap_idx[4]];
og_gad = aligned_masks[gad_idx[1]] .| aligned_masks[gad_idx[2]];
##
raphe_stim_exclude = L.exclude_from_masks([og_raphe], stim_exclude_masks[:,:,:,2])[1];
Gray.(stim_exclude_masks[:,:,10,2])

Gray.(og_raphe[:,:,10])
Gray.(raphe_stim_exclude[:,:,10])

## raphe atlas alignment looks suspicious...
img = RGB.(zeros(UInt8,size(og_raphe[:,:,10])...));
channelview(img)[[1],:,:] .= reshape(og_raphe[:,:,10],1,512,512)
# channelview(img)[[2,3],:,:] .= reshape(raphe_stim_exclude[:,:,10],1,512,512)
channelview(img)[[2,3],:,:] .= reshape(stim_exclude_masks[:,:,10,2],1,512,512)
img

## gad cluster 2
img = RGB.(zeros(UInt8,size(og_gad[:,:,1])...));
channelview(img)[[1],:,:] .= reshape(og_gad[:,:,9],1,512,512)
channelview(img)[[2,3],:,:] .= reshape(stim_exclude_masks[:,:,10,1],1,512,512)
img
##
# ex_mask = per_stim_region_masks[rap_idx[3],2] .& per_stim_region_masks[rap_idx[4],2]
ex_mask = per_stim_region_masks[rap_idx[3],2] .| per_stim_region_masks[rap_idx[4],2]
@show sum(ex_mask)
Gray.(ex_mask[:,:,10])
## debug on raphe 
regions_df = L.per_trial_regions_df_exclude_targets(
    tseries, window_len,
    stim_start_idx, stim_end_idx,
    # ss, se,
    trial_order, per_stim_region_masks[rap_idx[[3,4]],:], region_names[rap_idx[[3,4]]])  # , stim_names=nothing)

@show rap_names[3], sum(per_stim_region_masks[rap_idx[3],2])
filter(r->r.region == "Rhombencephalon - Raphe - Superior",
    regions_df)
## compare to no exclusion values...
no_exclude_masks =  repeat(aligned_masks[rap_idx[[3,4]]], 1, nStimuli)
regions_no_exclude_df = L.per_trial_regions_df_exclude_targets(
    tseries, window_len,
    stim_start_idx, stim_end_idx,
    # ss, se,
    trial_order, no_exclude_masks, region_names[rap_idx[[3,4]]])  # , stim_names=nothing)


##
@assert length(trial_order) == length(stim_start_idx)
# ss = stim_start_idx[1:length(trial_order)]
# se = stim_end_idx[1:length(trial_order)]
# @assert length(per_stim_region_masks) == length(region_names)
# TODO: exclude stim mask for appropriate stim
regions_df = L.per_trial_regions_df_exclude_targets(
    tseries, window_len,
    stim_start_idx, stim_end_idx,
    # ss, se,
    trial_order, per_stim_region_masks, region_names)  # , stim_names=nothing)
##
# Path to write out to (in the case that tyh5_path is defined, otherwise reading from tiff's)
out_dir = joinpath(fish_dir, basename(recording[:tyh5_path]))
out_dir = replace(out_dir, ".h5" => "")
# out_dir = joinpath(fish_dir, basename(recording_name), "raw")
println("Will create output_dir='$out_dir'")

mkpath(out_dir)
using Arrow
open(joinpath(out_dir, "regions_df.arrow"), "w") do io
    Arrow.write(io, regions_df)
end