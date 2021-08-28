using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
    Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
    NRRD, NIfTI
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
init_workers()
##
resource = Resources()
# TODO: also plot fish2..?
r = Recordings["2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"](
    ;resource...);
@pun (tseries, tseriesT, nStimuli, etl_vals, cells) = r;
## quick hack to fix bad recording params, should be fixed...
# r.nodes[:mm_warp_prefix] = "20210729T154219098Z"
# @pun (mm_warp_prefix, oir_dir) = r
# multimap_warpedname = glob_one_file("$mm_warp_prefix*SyN_Warped.nii.gz", oir_dir)
@pun oir_820_name = r
mm920_registered = niread(multimap_warpedname)
mm_transform_affine = glob_one_file("$mm_warp_prefix*GenericAffine.mat", oir_dir)
mm_transform_SyN = glob_one_file("$mm_warp_prefix*SyN_1Warp.nii.gz", oir_dir)
mm_transforms = [mm_transform_affine, mm_transform_SyN]
oir_820_file = joinpath(oir_dir, oir_820_name)
multimap_820 = L.read_olympus(oir_820_file);
mm820_registered = L.multimap_transforms(zseries, multimap_820,
    mm_transform_affine, mm_transform_SyN);


##
@pun (trial_average, vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, mm820_registered, mm920_registered, window_secs
) = r
##
@pun (stim_start_idx, stim_end_idx, vol_rate, zseries, plot_dir,
    tseriesZ, tseriesH, tseriesW, imaging2zseries_plane, target_size_px,
    recording_folder, fish_name, exp_name
) = r;
# (7056, 1, 10, 512, 512)
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)

# imap_cmap = transparent_cmap(plt.cm.Reds, max_alpha=1.0)
imap_cmap = transparent_cmap(plt.cm.Greens, max_alpha=1.0)
# bg_cmap = transparent_cmap(plt.cm.gray, max_alpha=0.5)
mm4_cmap = transparent_cmap(plt.cm.Reds, max_alpha=0.5)
mm1_cmap = transparent_cmap(plt.cm.Reds, max_alpha=0.5)
# mm4_cmap = transparent_cmap(plt.cm.Blues, max_alpha=0.5)
# mm1_cmap = transparent_cmap(plt.cm.Greens, max_alpha=0.5)
# TODO: plot I-map & multiMAP on same plot

max_frames = stim_start_idx[1]-1 
@warn "TODO: double check these"
window = minimum([Int(ceil(3*vol_rate)), max_frames])
# @assert (window <= post) & (window <= pre)
# @assert Z == 1
tseries_cmax = 2.5
tseries_cmin = -0.5

# tseries_cmax = 1.0
# tseries_cmax = 0.25 # fish2
tseries_cmax = 0.5 # fish1
tseries_cmin = 0.0

# mm4_cmax = quantile(mm820_registered[:,:,4,:][:], .999)
# mm4_cmax = quantile(ez_gamma(mm820_registered[:,:,4,:])[:], .99)
# mm4_cmin = otsu_threshold(ez_gamma(mm820_registered[:,:,4,:])[:])
mm4_cmin = otsu_threshold(mm820_registered[:,:,4,:][:]) * 2
mm4_cmax = quantile(mm820_registered[:,:,4,:][:], .99)
# mm4_cmin = 0
# mm4_cmax = 1

# mm4_cmax = 10
# mm1_cmax = quantile(mm820_registered[:,:,1,:][:], .999)
# mm4_cmin = 0.0
# mm1_cmax = quantile(ez_gamma(mm820_registered[:,:,1,:])[:], .99)
mm1_cmin = otsu_threshold(mm820_registered[:,:,1,:]) * 2
mm1_cmax = quantile(mm820_registered[:,:,1,:][:], .99)
# mm1_cmax = 1
# mm1_cmax = 10
# mm1_cmin = 0.3 # otherwise too much background

# tseries_cmax = 0.5
# tseries_cmin = -0.1
# tseries_cmax = 4
# tseries_cmin = -0.75
# tseries_cmin = -4
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=tseries_cmin,vcenter=0,vmax=cmax)

analysis_name = "lstm-multiMAP-zoom2x-xmlmatch"
# mm = 1
mm = 4

@warn "df_f denominator epsilon may have changed"
plus_minus = 100
yrange = (256-plus_minus, 256+plus_minus)
yrange = yrange .* 2
xrange = (225-plus_minus, 225+plus_minus)
xrange = xrange .* 2
##
for stimNum in 1:nStimuli
# for stimNum in [2]
    # f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    # TODO: should window be immediately after stim..? right now it's
    # ~2 seconds after stim...?
    # f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    # f0 = mean(trial_average[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    f0 = mean(trial_average[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    df = f - f0
    # originally no epsilon prior to July 2021
    # df_f = df./(f0 .+ 0.02)
    df_f = df./(f0 .+ 0.001)
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    if tseriesZ > 5
        figW,figH = (figB*tseriesZ/2, figB*2)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(2,Int(tseriesZ/2), figsize=(figW,figH))
    else
        figW,figH = (figB*tseriesZ, figB)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(1,tseriesZ, figsize=(figW,figH))
    end
    if tseriesZ==1
        axs = [axs]
    end
    # ax = axs[2]
    # ax = axs
    for z in 1:tseriesZ
        ax = axs[z]
        # global cim = ax.imshow(df_f[:,:,z], cmap="RdBu_r",
        #     norm=cnorm)
        match_z = imaging2zseries_plane[z]
        # fix scaling issue
        df_size = size(trial_average)[1:2]
        zseriesSize = size(zseries)[1:2]
        # influence map
        
        # global cim = ax.imshow(df_f[:,:,z], cmap=imap_cmap,
        #     clim=(tseries_cmin,tseries_cmax))
        
        # global cim = ax.imshow(imresize(df_f[:,:,z], zseriesSize), cmap=imap_cmap,
        #     clim=(tseries_cmin,tseries_cmax))
        # zseries
        # zs = ez_gamma(imresize(zseries[:,:,match_z], df_size))
        zs = ez_gamma(zseries[:,:,match_z])
        im = RGB.(imadjustintensity(zs))
        channelview(im)[[1,3],:,:] .= 0
        # ax.imshow(zs, cmap=bg_cmap)
        # multiMAP
        # ax.imshow(mm820_registered[:,:,4,match_z], cmap=mm4_cmap,
        #     clim=(mm4_cmin,mm4_cmax))
        # ax.imshow(mm820_registered[:,:,1,match_z], cmap=mm1_cmap,
        #     clim=(mm1_cmin,mm1_cmax))
        # ax.imshow(ez_gamma(imresize(mm820_registered[:,:,4,match_z], df_size)), cmap=mm4_cmap,
        #     clim=(mm4_cmin,mm4_cmax))
        # ax.imshow(ez_gamma(imresize(mm820_registered[:,:,1,match_z], df_size)), cmap=mm1_cmap,
        #     clim=(mm1_cmin,mm1_cmax))
        # ax.imshow(ez_gamma(mm820_registered[:,:,4,match_z]), cmap=mm4_cmap,
        #     clim=(mm4_cmin,mm4_cmax))
        if mm == 4
            # ax.imshow(mm820_registered[:,:,4,match_z], cmap=mm4_cmap, # use this one
            #     clim=(mm4_cmin,mm4_cmax)) # use this one
            channelview(im)[1,:,:] .= imadjustintensity(mm820_registered[:,:,4,match_z],
                (mm4_cmin,mm4_cmax))
            ax.imshow(im)
        end
        # ax.imshow(ez_gamma(mm820_registered[:,:,1,match_z]), cmap=mm1_cmap,
        #     clim=(mm1_cmin,mm1_cmax))
        if mm==1
            # ax.imshow(mm820_registered[:,:,1,match_z], cmap=mm1_cmap, # use this one
            #     clim=(mm1_cmin,mm1_cmax)) 
        end
        ax.set_axis_off()
        # TODO: fix error on next line
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_title("$(Int(round(etl_vals[z],digits=0)))μm")
        # this will make extra circles (1 extra per repetition...)
        for (x,y,targetZ) in eachrow(unique(cells[cells.stimNum .== stimNum,[:x,:y,:z]]))
            if z == targetZ
                circle = matplotlib.patches.Circle((x,y), target_size_px, color="k",
                    fill=false, lw=0.4, alpha=0.3)
                ax.add_patch(circle)
            end
        end
        # previous (but not this stim) targets
        # for (x,y,z) in eachrow(unique(cells[cells.stimNum .!= stimNum,[:x,:y,:z]]))
        #     @assert z == 1
        #     circle = matplotlib.patches.Circle((x,y), targetSizePx, color="k",
        #         fill=false, lw=0.5, alpha=0.5)
        #     ax.add_patch(circle)
        # end
    end

    # axs[1].imshow(stim_masks[:, :,z,stimNum], cmap="gray")

    # cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
    # plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
    #     norm=cnorm)
    # may need to adjust if colorbar is cutoff
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.15, 0.0075, 0.7])
    # cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
    cbar = fig.colorbar(cim, cax=cbar_ax)
    path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_stim$(stimNum)_mm$mm")
    @show path*".svg"
    fig.savefig(path*".svg", dpi=1200)
    fig.savefig(path*".png", dpi=1200)
end