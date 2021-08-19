ENV["DISPLAY"] = "localhost:10"
##
using AlgebraOfGraphics, CairoMakie
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays
import Lensman: @pun, @assign
using ImageView
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
@pun h2b_zbrain = resources
imshow(h2b_zbrain)


## TODO: takes 20 seconds due to precompilation

r = Recordings[
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    ](;resources...);

##
@pun (region_masks_h5, imaging2zseries_plane,
    zbrain_masks, zbrain_mask_names, plot_dir, 
    tseries, trial_average, window_len, trial_order,
    # stim_start_idx, stim_end_idx,
    tyh5_path, tseries_dset, 
    stimStartIdx, stimEndIdx,
    tseriesH, tseriesW, tseriesZ
) = r;
stim_start_idx, stim_end_idx = stimStartIdx, stimEndIdx
##
# @pun (nCells, cell_centers, cells_mask, iscell) = r
# @pun suite2p_dir = r

##
imap = influence_map(trial_average, window_len);
##
function read_first_mask(mask_name)
    name = zbrain_mask_names[occursin.(mask_name, zbrain_mask_names)][1]
    name, read_registered_mask(region_masks_h5,
        name)[:,:,imaging2zseries_plane];
end

function read_matching_mask(mask_name)
    names = zbrain_mask_names[occursin.(mask_name, zbrain_mask_names)]
    ret = nothing
    for name in names
        new = read_registered_mask(region_masks_h5, name)[:,:,imaging2zseries_plane]
        if isnothing(ret)
            ret = new
        else
            ret = ret .+ new
        end
    end
    names, ret
end
tectum_name, tectum = read_first_mask("Tectum")
tectum_name
cerebellum_names, cerebellum = read_matching_mask("Cerebe")
cerebellum_names
medulla_names, medulla = read_matching_mask("Medulla")
medulla_names
##
raphe_names, raphe = read_matching_mask("Raphe")
raphe_names

##
habenula_name, habenula = read_first_mask("Habenula")
habenula_name
rhom1_name, rhom1 = read_first_mask("Rhombomere 1")
rhom1_name

stim_names = ["left habenula", "right habenula", "control"]
region_names = ["habenula", "optic tectum", "rhombomere 1", "raphe", "cerebellum"]
masks = [habenula, tectum, rhom1, raphe, cerebellum];
masks = map(m->imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);
# regions = [habenula_name, tectum_name, rhom1_name]
# stim_names = ["left habenula", "right habenula", "control"]
# pretty_region_names = ["habenula", "optic tectum", "rhombomere 1"]

##
"Return DataFrame of df/f values"
function imap_DF_by_region(imap, region_masks_h5::HDF5.File, regions, imaging2zseries_plane,
    stim_names=nothing, pretty_region_names=nothing
)
    nStimuli = size(imap,4)
    if isnothing(stim_names)
        stim_names = string.(collect(1:nStimuli))
    end
    pretty_region_names = isnothing(pretty_region_names) ? regions : pretty_region_names

    H, W, Z = size(imap)[1:3]
        
    df = DataFrame()
    for (pretty, region) in zip(pretty_region_names, regions)
        mask = read_registered_mask(region_masks_h5, region)[:,:,imaging2zseries_plane]
        for s in 1:nStimuli
            voxel_responses = selectdim(imap,4,s)[mask]
            voxel_responses = mean(voxel_responses)
            N = length(voxel_responses)
            df = vcat(df, DataFrame("Δf/f"=>voxel_responses,
                # region=repeat(CategoricalArray([region]),N),
                "region"=>pretty,
                # stim=repeat(CategoricalArray([stim_names[s]]),N)))
                "stim"=>stim_names[s]))
        end
    end
    df
end

function imap_DF_by_region(imap, masks, regions, imaging2zseries_plane,
    H, W, Z, stim_names=nothing
)
    nStimuli = size(imap,4)
    if isnothing(stim_names)
        stim_names = string.(collect(1:nStimuli))
    end

    H, W, Z = size(imap)[1:3]
        
    df = DataFrame()
    for (region, mask) in zip(regions, masks)
        for s in 1:nStimuli
            voxel_responses = selectdim(imap,4,s)[mask]
            N = length(voxel_responses)
            df = vcat(df, DataFrame("Δf/f"=>voxel_responses,
                "region"=>region,
                "stim"=>stim_names[s]))
        end
    end
    df
end

@everywhere include("/home/tyler/code/lensman/src/tidying.jl")
@everywhere using DataFrames, Statistics
regions_df = per_trial_regions_df(tyh5_path, tseries_dset, window_len,
    stim_start_idx, stim_end_idx,
    # stimStartIdx, stimEndIdx,
    trial_order, masks, region_names, stim_names)

##
# p1 = Data(regions_df) * visual(Violin) *
p1 = Data(regions_df) * visual(BoxPlot, show_notch=true) *
    mapping(:region, "Δf/f", color=:stim, dodge=:stim) |> aog.draw

ppath = joinpath(plot_dir, "region_df_f_by-trial.png")
save(ppath, p1)
##

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







##
Gray.(tectum[:,:,5])
##
Gray.(imadjustintensity(mean(trial_average, dims=[4,5])[:,:,5,1,1]))

##

Gray.(habenula[:,:,6])
##

sum(tectum)

Gray.(tectum[:,:,10])

sum(tectum, dims=[1,2])


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