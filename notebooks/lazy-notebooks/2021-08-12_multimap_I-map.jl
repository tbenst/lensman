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
##
# TODO: also plot fish2..?
r = Recordings[
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    ](
    ;resource...);
##
# analysis_name = "lstm-multiMAP-zoom3x-2021-08-28"
analysis_name = "multiMAP-2021-08-28"
# analysis_name = "lstm-multiMAP-zoom1x-2021-08-28"

##
@pun (tseries, tseriesT, nStimuli, etl_vals, cells,
    mm920_registered, mm820_registered
) = r;
@pun (trial_average, vol_rate, window_len, recording_folder, fish_name,
    exp_name, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, mm820_registered, mm920_registered, window_secs
) = r
@pun (stim_start_idx, stim_end_idx, vol_rate, zseries, plot_dir,
    tseriesZ, tseriesH, tseriesW, imaging2zseries_plane, target_size_px,
    recording_folder, fish_name, exp_name
) = r;
##
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

# tseries_cmax = 
tseries_cmax = 0.25 # fish2
tseries_cmax = quantile(trial_average[:], .99)
# tseries_cmax = 0.75 # fish1
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
mm1_cmax = quantile(mm820_registered[:,:,1,:][:], .999)
mm2_cmin = otsu_threshold(mm820_registered[:,:,2,:]) * 2
mm2_cmax = quantile(mm820_registered[:,:,2,:][:], .999)
mm3_cmin = otsu_threshold(mm820_registered[:,:,3,:]) * 2
mm3_cmax = quantile(mm820_registered[:,:,3,:][:], .999)
# mm1_cmax = maximum(mm820_registered[:,:,1,:][:])
# mm1_cmax = 1
# mm1_cmax = 10
# mm1_cmin = 0.3 # otherwise too much background

# tseries_cmax = 0.5
# tseries_cmin = -0.1
# tseries_cmax = 4
# tseries_cmin = -0.75
# tseries_cmin = -4
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=tseries_cmin,vcenter=0,vmax=cmax)
@warn "df_f denominator epsilon may have changed"
plus_minus = 100
##
## plot one plane for quick iteration
yrange = (256-plus_minus, 256+plus_minus)
# yrange = yrange .* 2
xrange = (225-plus_minus, 225+plus_minus)
# xrange = (256-plus_minus, 256+plus_minus)
# xrange = xrange .* 2
mm=4
stimNum = 1


f = mean(trial_average[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
f0 = mean(trial_average[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
df = f - f0
# originally no epsilon prior to July 2021
# df_f = df./(f0 .+ 0.02)
df_f = df./(f0 .+ 0.001)
figW,figH = (figB*1, figB)
global fig = plt.figure(figsize=(figW,figH))
fig, ax = plt.subplots(1,1, figsize=(figW,figH))
z = 5


match_z = imaging2zseries_plane[z]
# fix scaling issue
df_size = size(trial_average)[1:2]
zseriesSize = size(zseries)[1:2]

# influence map


# im = imadjustintensity(df_f[:,:,z], (tseries_cmin,tseries_cmax))
im = df_f[:,:,z]
im = opening_nonbool(im)
cmin = quantile(im[:], 0.01)
cmax = quantile(im[:], 0.99)
im = imadjustintensity(im, (cmin, cmax))
# im = mapwindow(median!, im, (3,3))
im = RGB.(im)
channelview(im)[[1,3],:,:] .= 0

# multimap
if mm==4
    cmin, cmax = mm4_cmin, mm4_cmax
elseif mm==1
    cmin, cmax = mm1_cmin, mm1_cmax
elseif mm==2
    cmin, cmax = mm2_cmin, mm2_cmax
end
im820 = imresize(mm820_registered[:,:,mm,match_z],
        df_size)
im820 = imadjustintensity(opening_median(im820))
channelview(im)[1,:,:] .= im820
# channelview(im)[1,:,:] .= 0




py_im = permutedims(channelview(im),[2,3,1])
cim = ax.imshow(py_im)
ax.set_axis_off()
ax.set_xlim(xrange)
ax.set_ylim(yrange)

fig

## check registration
z = 1
# mm=4
mm=3
match_z = imaging2zseries_plane[z]
fig, ax = plt.subplots(1,1, figsize=(figW,figH))
zplane = imresize(Array(zseries[:,:,imaging2zseries_plane][:,:,z]), df_size)
im = imadjustintensity(zplane)
# im = imadjustintensity(zplane, (tseries_cmin,tseries_cmax))
im = opening_nonbool(im)
im = mapwindow(median!, im, (3,3))
# im = mapwindow(median!, im, (3,3))


# im = imadjustintensity(imresize(mm920_registered[:,:,match_z],
#     df_size))


im = RGB.(im)
im = adjust_histogram(im, GammaCorrection(0.5))
channelview(im)[[1,3],:,:] .= 0
channelview(im) .= 0

# multimap
if mm==4
    cmin, cmax = mm4_cmin, mm4_cmax
elseif mm==1
    cmin, cmax = mm1_cmin, mm1_cmax
elseif mm==2
    cmin, cmax = mm2_cmin, mm2_cmax
elseif mm==3
    cmin, cmax = mm3_cmin, mm3_cmax
end
# im820 = imadjustintensity(imresize(mm820_registered[:,:,mm,match_z],
#     df_size))
# im820 = imadjustintensity(imresize(mm820_registered[:,:,mm,match_z],
#     df_size), (cmin, cmax))
    # df_size), (0,0))
    # df_size), (cmin, cmin))
# im820 = imadjustintensity(imresize(mm920_registered[:,:,match_z],
#     df_size))
# im820 = opening_median(im820)

# im820 = imresize(mm820_registered[:,:,mm,match_z],
#         df_size)
im820 = opening_median(im820)
# cmin, cmax = otsu_threshold(im820[:])*2, quantile(im820[:], 0.99)
# im820 = adjust_histogram(im820, GammaCorrection(0.5))
cmin, cmax = minimum(im820[:]), minimum(im820[:])
im820 = imadjustintensity(im820, (cmin, cmax))



# im820 .= 0
channelview(im)[1,:,:] .= im820

py_im = permutedims(channelview(im),[2,3,1])
cim = ax.imshow(py_im)
ax.set_axis_off()
ax.set_xlim(xrange)
ax.set_ylim(yrange)

fig

##
# mm=4
mm=1
for stimNum in nStimuli:-1:1
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
        match_z = imaging2zseries_plane[z]
        # fix scaling issue
        df_size = size(trial_average)[1:2]
        zseriesSize = size(zseries)[1:2]
        
        # influence map
        im = df_f[:,:,z]
        im = opening_median(im)
        cmin, cmax = quantile(im[:], 0.01), quantile(im[:], 0.99)
        im = imadjustintensity(im, (cmin, cmax))
        im = RGB.(im)
        channelview(im)[[1,3],:,:] .= 0

        # multimap
        if mm==4
            cmin, cmax = mm4_cmin, mm4_cmax
        elseif mm==1
            cmin, cmax = mm1_cmin, mm1_cmax
        end
        im820 = imresize(mm820_registered[:,:,mm,match_z],
            df_size)
        cmin = otsu_threshold(im820[:])*2
        cmax = quantile(im820[:], 0.99)
        im820 = imadjustintensity(im820, (cmin, cmax))
        im820 = opening_median(im820)

        channelview(im)[1,:,:] .= im820

        py_im = permutedims(channelview(im),[2,3,1])
        cim = ax.imshow(py_im)

        ax.set_axis_off()
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_title("$(Int(round(etl_vals[z],digits=0)))μm")
    end

    path = joinpath(plot_dir,"$(recording_folder)_$(fish_name)_$(exp_name)_$(analysis_name)_stim$(stimNum)_mm$mm")
    @show path*".svg"
    fig.savefig(path*".svg", dpi=1200)
    fig.savefig(path*".png", dpi=1200)
end
fig