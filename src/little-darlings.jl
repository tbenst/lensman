# code that is not used but I can't bear parting with

function adj_gamma_rgb(img; gamma=0.5)
    colorview(RGB, adjust_histogram(channelview(img), GammaCorrection(gamma)))
end

"combine three seperate 1-channel RGB into a 3-channel RGB."
function singlechanrgb2rgb(red, green, blue)
    new = copy(red)
    greenview = selectdim(channelview(new),1,2)
    greenview[:] = selectdim(channelview(green),1,2)
    blueview = selectdim(channelview(new),1,3)
    blueview[:] = selectdim(channelview(blue),1,3)
    new
end

"layer ontop of base."
function layer_imgs(base, layer; ratio=0.8)
    new = base .* (1-ratio)
    new += layer .* ratio
    new
end

"""Compute and assign local variables.
```
@experiment (etlVals, Z, H, W) = "2021-08-03/fish1/TSeries-1"
```

will expand to:
```
etlVals, Z, H, W = compute("2021-08-03/fish1/TSeries-1",
    ["etlVals", "Z", "H", "W"])
```
"""
macro experiment(ex)
    @assert ex.head == :(=) "no `=` found in expression."
    vars = ex.args[1] # :((etlVals, Z, H, W))
    exp_name = ex.args[2] # "2021-08-03/fish1/TSeries-1"
    vars_str = map(string, vars.args) # ["etlVals", "Z", "H", "W"]
    esc(:($vars = compute($exp_name, $vars_str)))
    # ex



## code exploring NWB .nwb file for suite2p output
using HDF5
# h5 = h5open("/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123.ty.h5", "r")
ophys = h5open("$suite2p_dir/ophys.nwb", "r")
# 2 x nNeurons
size(ophys["/processing/ophys/ImageSegmentation/PlaneSegmentation/iscell"])
# (2, 11670)

size(ophys["/processing/ophys/ImageSegmentation/PlaneSegmentation/voxel_mask"])
# (556770,) => (x = 0x00000000, y = 0x000000b4, z = 0x00000000, weight = 0.0014966136f0)


ophys["/processing/ophys/ImageSegmentation/"]
# nNeurons
size(ophys["/processing/ophys/ImageSegmentation/PlaneSegmentation/voxel_mask_index"])
# (11670,) => uint, all unique numbers from (421, 556770)

# nTime x nNeurons
size(ophys["/processing/ophys/Fluorescence/Fluorescence/data"])
# (11444, 11670)

# worthless index list..?
size(ophys["/processing/ophys/Fluorescence/Fluorescence/rois"])


##
    
## pattern for extending DAG (ugly! ergonomics should be improved..macro?)
global df_f_per_trial_dataframe
if ~in(:df_f_per_trial_dataframe, keys(r.thunks))
    @lazy begin
        window_len = ((vol_rate)->Int(floor(5 * vol_rate)) - 1)(vol_rate)
        df_f_per_trial_dataframe = get_df_f_per_trial_dataframe(
            df_f_per_voxel_per_trial, trial_order)
    end
    @assign r.thunks = df_f_per_trial_dataframe
end
@pun df_f_per_trial_dataframe = r;
