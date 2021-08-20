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



############## dead code from 2021-08-18_hemisphere_zbrain.jl

# stim_names = nothing
smallmasks= 
df = imap_DF_by_region(imap, masks, region_names, imaging2zseries_plane,
    512, 512, 151, stim_names)
# df = imap_DF_by_region(imap, region_masks_h5, regions, imaging2zseries_plane,
#     stim_names, pretty_region_names)
##
# plot(df, x="df_f", color="region", Geom.density, Scale.y_log10)
# Data(df) * visual(Violin) * mapping(:region, :df_f) |> aog.draw
# categorical: 37s
# strings: 18s
# one char string: 1.8s
# one char CategoricalArray, string: 6.8s
# two col, long CategoricalArray: 35s
# p0 = Data(df[1:100,:]) * visual(Violin) *
    # mapping(:region, :df_f, color=:stim, dodge=:stim)
p1 =  Data(df) * visual(Violin) *
    mapping(:region, "Δf/f", color=:stim, dodge=:stim)
    
# aog.draw(p0, axis=(xscale=log,))
# aog.draw(p0, axis=(xscale=log,))
# aog.draw(p1, axis=(xscale=log,))
d1 = aog.draw(p1)
ppath = joinpath(plot_dir, "voxel_df_f_by-region.png")
save(ppath, d1)
# aog.draw(p1; axis=(xscale=log,))

##
# Gadfly.plot(df, x="df_f", color="region", Geom.density, Scale.y_log10)
# Gadfly.plot(df, x="df_f", color="region", Geom.violin, Scale.y_log10)
Gadfly.plot(df, x=:stim, y=:df_f, color=:region,
    Geom.violin)






## not finished
figB = 1.6
figW,figH = (figB*Z, figB)
fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(1,Z, figsize=(figW,figH))
imap_cmap = transparent_cmap(plt.cm.Greens, max_alpha=1.0)
for z in 1:Z
    ax = axs[z]
    match_z = imaging2zseries_plane[z]
    df_size = size(trial_average)[1:2]
    zseriesSize = size(zseries)[1:2]
    
    global cim = ax.imshow(imresize(df_f[:,:,z], zseriesSize), cmap=imap_cmap,
        clim=(tseries_cmin,tseries_cmax))
    zs = ez_gamma(zseries[:,:,match_z])
    # ax.imshow(zs, cmap=bg_cmap)
    # ax.imshow(mm820_registered[:,:,4,match_z], cmap=mm4_cmap,
    #     clim=(mm4_cmin,mm4_cmax))
    # ax.imshow(mm820_registered[:,:,1,match_z], cmap=mm1_cmap,
    #     clim=(mm1_cmin,mm1_cmax))
    ax.set_axis_off()
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.set_title("$(Int(round(etlVals[z],digits=0)))μm")
    # this will make extra circles (1 extra per repetition...)
    for (x,y,targetZ) in eachrow(unique(cells[cells.stimNum .== stimNum,[:x,:y,:z]]))
        if z == targetZ
            circle = matplotlib.patches.Circle((x,y), targetSizePx, color="k",
                fill=false, lw=0.4, alpha=0.3)
            ax.add_patch(circle)
        end
    end
end

# may need to adjust if colorbar is cutoff
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.91, 0.15, 0.0075, 0.7])
cbar = fig.colorbar(cim, cax=cbar_ax)
# path = joinpath(plotDir,"$(recording_folder)_$(fish_name)_$(expName)_$(analysis_name)_stim$stimNum")
# @show path*".svg"
# fig.savefig(path*".svg", dpi=1200)
# fig.savefig(path*".png", dpi=1200)
fig

#################