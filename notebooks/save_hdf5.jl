# need X11 forwarding for ssh, `ssh -X osprey`, then `echo $DISPLAY`, and fill below
ENV["DISPLAY"] = "localhost:10.0"
using Lensman

## Define experiment paths

# TODO(allan.raventos): verify that writing without stim masks still works
# fish_dir = "/scratch/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2"
# exp_name = "TSeries-resting-pre-043"
# output_path = joinpath(fish_dir, exp_name, "$(exp_name).ty.h5")

fish_dir = "/scratch/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2"
exp_name = "TSeries-lrhab_raphe_40trial-045"
tif_dir = joinpath(fish_dir, exp_name)

## Load tseries
tseries = loadTseries(tif_dir)
(H, W, Z, T) = size(tseries)

## Get stim_masks and stim index for each timestep
# slm_dir = "/mnt/deissero/users/tyler/b115/SLM_files"
slm_dir = "/scratch/allan/b115"
z_offset = 0.
(stim_masks, stim_used_at_each_timestep) = get_slm_stim_masks(tif_dir, slm_dir, tseries, z_offset)
size(stim_masks)

## Write out experiment
output_path = "$(exp_name).ty.h5"
write_experiment_to_tyh5(tseries, stim_masks, stim_used_at_each_timestep, output_path)
