ENV["DISPLAY"] = "localhost:11"
##
using AlgebraOfGraphics, CairoMakie # should import first
using Lensman, Images, HDF5, Statistics, Colors, ColorVectorSpace,
    ColorTypes
import Lensman: @pun, @assign
L = Lensman
using Thunks, ImageView

##
resources = Resources();
##
@pun (zbrain_dir, zbrain_masks) = resources
options = Dict(
)
recording = Recordings[
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    ](;options...,
    resources...);

##
@pun (fish_dir, zseries, zbrain_mask_names, imaging2zseries_plane,
    plot_dir
) = recording;

##
region_mask_path = "$fish_dir/h5_output/region_masks.h5"
region_masks_h5 = h5open(region_mask_path, "r", swmr=true);
##
tectum_name, tectum = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Tectum")
cerebellum_names, cerebellum = L.read_matching_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Cerebe")
raphe_names, raphe = L.read_matching_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Raphe")
habenula_name, habenula = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Habenula")
rhom1_name, rhom1 = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombomere 1");
rhom2_name, rhom2 = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombomere 2");
rhom3_name, rhom3 = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombomere 3");
rhom4_name, rhom4 = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombomere 4");
rhom5_name, rhom5 = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombomere 5");
ipn_name, ipn = L.read_first_mask(region_masks_h5, zbrain_mask_names,
    imaging2zseries_plane,"Rhombencephalon - Interpeduncular Nucleus");
##
v = quantile(zseries[:,:,imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:,:,imaging2zseries_plane], (0,v))
thresh = otsu_threshold(green[:]) * 0.5
for c in findall(green .<= thresh)
    green[c] = 0
end
# green = adjust_histogram(green,GammaCorrection(0.5))
regions = [
    tectum,
    # cerebellum,
    # raphe,
    habenula,
    # ipn
    rhom1,
    rhom2,
    rhom3,
    rhom4,
    rhom5,
]
# we want to distinguish from any value of green
seeds = vcat([RGB(0,x,0) for x in 0:0.11:1]..., RGB(1,1,1), RGB(0,0,0))
colors = distinguishable_colors(length(regions), seeds, dropseed=true)


im2 = RGB.(green)
channelview(im2)[[1,3],:,:,:] .= 0
for (r, reg) in enumerate(regions)
    region = reg .> 0
    color_mask = RGB{N0f16}.(zeros(size(im2)...))
    # color_mask .= RGB{N0f16}(1,1,1)
    color_mask[region] .= colors[r] * 0.8
    im2 = mapc.((a,b)->maximum([a,b]), im2, color_mask)
end
# channelview(im2)[1,:,:,:] .= habenula .* 0.8
# channelview(im2)[1,:,:,:] .+= rhom1 .* 0.8
# channelview(im2)[3,:,:,:] .= tectum .* 0.8
##
# h_profile = sum(im2[:,:,5], dims=[2])[:,1]
save(joinpath(plot_dir, "hab_plane_region_annotations.png"), im2[:,:,4])
##
ipn_all = read_registered_mask(region_masks_h5,
        "Rhombencephalon - Interpeduncular Nucleus");
raphe_all = read_registered_mask(region_masks_h5,
        "Rhombencephalon - Raphe - Superior");

im1 = RGB.(ipn_all)
channelview(im1)[[1,3],:,:,:] .= 0
channelview(im1)[1,:,:,:] .= raphe_all
imshow(im1)
##
# todo: automatically crop?
vals = sum(ipn_all, dims=[1,2])[1,1,:]
mapping(1:length(vals), vals) |> draw