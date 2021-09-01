##
resources = Resources(zbrain_dir="/data/b115/atlas/zbrain");

##
@pun (zbrain_dir, zbrain_masks) = resources
options = Dict(
    :tseries_dset => nothing, :rel_analysis_dir => "plots"  # use raw data
)
recording = Recordings[
        # "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    ](;options...,
    resources...);

##
# @pun (region_masks_h5, imaging2zseries_plane,
#     zbrain_masks, zbrain_mask_names, plot_dir, 
#     tseries, trial_average, window_len, trial_order,
#     stim_start_idx, stim_end_idx, target_groups,
#     tyh5_path, tseries_dset, zbrain_transforms,
#     tseriesH, tseriesW, tseriesZ, zseries, rostral, dorsal,
#     recording_folder, fish_name,
#     exp_name
# ) = recording;
@pun (imaging2zseries_plane) = recording;

##
# file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/zseries.nrrd"
# zseries = load(file_path)

# Load atlas
file_path = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/zbrain_h2b.nrrd"
fish_atlas = load(file_path)

# Load baseline
using NIfTI
file_path = "/data/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/20210726T095924956Z_Warped.nii.gz"
baseline = niread(file_path)

@pun 
v = quantile(zseries[:, :, imaging2zseries_plane][:], 0.99)
green = imadjustintensity(zseries[:, :, imaging2zseries_plane], (0, v))
im2 = RGB.(green)
channelview(im2)[1,:,:,:] .= fish_atlas .* 0.8
imshow(im2)
