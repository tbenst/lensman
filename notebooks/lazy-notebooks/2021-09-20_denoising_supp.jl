ENV["DISPLAY"] = "localhost:11"
##
# using ImageView
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Unitful
import PyPlot
import Lensman: @pun, @assign
plt = PyPlot
matplotlib = plt.matplotlib

SMALL_SIZE = 7
MEDIUM_SIZE = 8
BIGGER_SIZE = 9
plt.rc("font", size=SMALL_SIZE)          # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

np = pyimport("numpy")
matscale = pyimport("matplotlib_scalebar.scalebar")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman

matplotlib.rcParams["font.sans-serif"] = ["Arial", "sans-serif"]

init_workers()
resources = Resources();
##
r_lstm = Recordings[
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...,
    # trial_average_path="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_lstm_avgStim.h5",
    # trial_average_dset="/lstm"
);
@pun (plot_dir, trial_average_path, trial_average_dset, tseries_units,
    window_len, region_masks_h5, zbrain_mask_names, imaging2zseries_plane,
    tseriesH, tseriesW, tseriesZ
) = r_lstm
# h5 = h5open(trial_average_path)
size(h5read(trial_average_path,trial_average_dset))
# size(h5read(trial_average_path,"lstm"))
##
r_kalman = Recordings[
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...,
    tseries_read_strategy = :lazy_hwzt,
    tyh5_path="/data/dlab/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5",
    tseries_dset=nothing
);
feng_lstm_path = "/data/dlab/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_fenglstm.ty.h5"
r_feng = Recordings[
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...,
    tseries_read_strategy = :lazy_tyh5,
    tyh5_path=feng_lstm_path,
    tseries_dset="/imaging/original_LSTM_with_trainable_initial_states-2021-07-14",
    trial_average_path = "/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_fenglstm_avgStim.h5"
);
r_raw = Recordings[
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...,
    tseries_read_strategy = :lazy_tiff,
);
r_feng[:trial_average_path]
##
projective_field_lstm = r_lstm[:projective_field];
projective_field_raw = r_raw[:projective_field];
projective_field_kalman = r_kalman[:projective_field];
println("feng1")
projective_field_feng = r_feng[:projective_field];
##
trial_average_lstm = r_lstm[:trial_average];
trial_average_raw = r_raw[:trial_average];
trial_average_kalman = r_kalman[:trial_average];
trial_average_feng = r_feng[:trial_average] .* 8192;

## get outlines
REGION_LIST = [
    "Rhombencephalon -"
    "Mesencephalon -"
    "Diencephalon -"
    "Telencephalon -"
]
region_masks = [L.read_first_mask(region_masks_h5, zbrain_mask_names,
imaging2zseries_plane, region; outline=false) for region in REGION_LIST];

region_outlines = [dilate(L.mask2outline(L.unsalt_n_pepper(regm[2]; felz_min=5000)),[1,2])
    for regm in region_masks];
    # for regm in region_masks[[1]]];
outlines = cat(region_outlines...,dims=4);
outlines = maximum(outlines,dims=4)[:,:,:,1];
small_outlines = imresize(outlines, (tseriesH, tseriesW, tseriesZ));
Gray.(outlines[:,:,3])


##
if false
tseries_feng = r_feng[:tseries];
tseries_raw = r_raw[:tseries];
tseries_kalman = r_kalman[:tseries];
# trial_average_feng = r_feng[:trial_average];
@show minimum(trial_average_feng), mean(trial_average_feng), maximum(trial_average_feng)
@show minimum(trial_average_raw), mean(trial_average_raw), maximum(trial_average_raw)

## variance
# http://www.csun.edu/~mr31841/documents/Varitionandpredictionintervals.pdf
avg_raw = Float64.(trial_average_raw)
avg_kalman = Float64.(trial_average_kalman)
avg_feng = Float64.(trial_average_feng)
# raw_tv = sum(var(avg_raw,dims=5)[:,:,:,:,1])
raw_tv = sum((avg_raw .- mean(avg_raw,dims=5)) .^ 2)
kalman_explained = sum((avg_kalman .- mean(avg_raw,dims=5)) .^ 2)
kalman_unexplained = sum((avg_kalman .- avg_raw) .^ 2)
feng_explained = sum((avg_feng .- mean(avg_raw,dims=5)) .^ 2)
feng_unexplained = sum((avg_feng .- avg_raw) .^ 2)
@show raw_tv, kalman_explained, kalman_unexplained, feng_explained, feng_unexplained
kalman_unexplained, feng_unexplained

norm(avg_feng-avg_raw)
norm(avg_kalman-avg_raw)
##
im1_feng = permutedims(tseries_feng[:,:,4,100], (2,1))
im2_feng = permutedims(tseries_feng[:,:,4,500], (2,1))
im3_feng = permutedims(tseries_feng[:,:,4,1000], (2,1))
im_feng = hcat(im1_feng,im2_feng,im3_feng) .* 8192
im1_raw = tseries_raw[:,:,4,100]
im2_raw = tseries_raw[:,:,4,500]
im3_raw = tseries_raw[:,:,4,1000]
im_raw = hcat(im1_raw,im2_raw,im3_raw)
im1_kalman = tseries_kalman[:,:,4,100]
im2_kalman = tseries_kalman[:,:,4,500]
im3_kalman = tseries_kalman[:,:,4,1000]
im_kalman = hcat(im1_kalman,im2_kalman,im3_kalman)
imgrid = vcat(im_raw,im_kalman, im_feng);
cmax = quantile(imgrid[:],0.95)
@show norm(im1_feng .- im1_raw), norm(im1_kalman .- im1_raw)
diff_kalman = (im2_kalman .- im1_kalman);
diff_feng = (Float64.(im2_feng) .- Float64.(im1_feng));
diff_raw = (im2_raw .- im1_raw);
@show norm(diff_feng .- diff_raw), norm(diff_kalman .- diff_raw)
Gray.(imadjustintensity(imgrid, (0,cmax)))
##
Gray.(imadjustintensity(hcat(im1,im2,im3)))

# tseries_feng
end # if false

## compare projective fields including feng
figB = Float64(uconvert(u"inch", 183u"mm") / u"inch")
figW = figB
figH = figB / 2
fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(2,4, figsize=(figW,figH))
z = 4
s = 16
cmax = 1.0
cmin = 0
# cmin = -0.5
cmap = "viridis"
# cmap="RdBu_r"
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
df_f = projective_field_raw[:,:,z,16]
df_f_nosalt = L.unsalt_n_pepper(df_f)
cim = axs[1].imshow(df_f, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[1,1].set_title("Raw")
cim = axs[2,1].imshow(df_f_nosalt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[2,1].set_title("Raw - felzenszwalb")

df_f = projective_field_kalman[:,:,z,16]
df_f_nosalt = L.unsalt_n_pepper(df_f)
axs[1,2].imshow(df_f, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[1,2].set_title("Kalman")
axs[2,2].imshow(df_f_nosalt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[2,2].set_title("Kalman - felzenszwalb")


df_f = projective_field_feng[:,:,z,16];
df_f_nosalt = L.unsalt_n_pepper(df_f)
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
cim = axs[1,4].imshow(df_f, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[1,4].set_title("LSTM per pixel")
axs[2,4].imshow(df_f_nosalt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[2,4].set_title("LSTM per pixel - felzenszwalb")


cmax = 0.02
cmin = -0.02


# cmax = 0.03
# cmin = -0.05

axs[1,3].imshow(projective_field_lstm[:,:,z,16], cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
df_f = projective_field_lstm[:,:,z,16]
df_f_nosalt = L.unsalt_n_pepper(df_f)
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
axs[1,3].imshow(df_f, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[1,3].set_title("shared LSTM w/ trainable init")
axs[2,3].imshow(df_f_nosalt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
    # norm=cnorm, interpolation=nothing)
axs[2,3].set_title("shared LSTM - felzenszwalb")

for ax in axs
    ax.set_axis_off()
end
plt.tight_layout()

analysis_name = "denoising_comparison"
plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_df_f")
fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=300)
@show plotpath*".svg"

fig
## Spatial denoising FINAL FIG
figB = Float64(uconvert(u"inch", 183u"mm") / u"inch")
figW = figB
figH = figB / 3
fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(1,4, figsize=(figW,figH))
z = 3
s = 16
cmax = 1.0
cmin = 0
# cmin = -0.5
# cmap = "viridis"
cmap = "inferno"
# df_f = projective_field_raw[:,:,z,16]
# df_f = projective_field_lstm[:,:,z,16]
# df_f = projective_field_kalman[:,:,z,16]

window = 9
f = mean(trial_average_lstm[:,:,z,16,end-window+1:end],dims=3)[:,:,1]
f0 = mean(trial_average_lstm[:,:,z,16,1:window],dims=3)[:,:,1]
df = f - f0
df_f = df ./ f0


# slice = 256-50:256+50, 300:400
# slice = 256-50:256+50, 256-50:256+50
slice = 1:512,1:512
df_f = df_f[slice...]

felz_k=50
# felz_k=100
felz_min=10
df_f_nosalt = L.unsalt_n_pepper(df_f, felz_k, felz_min)
opened = opening_nonbool(df_f)
med_filt = opening_median(df_f)
segments = felzenszwalb(med_filt, felz_k, felz_min)
salt = copy(opened)
salt[segments.image_indexmap .!= 1] .= 0
# salt[segments.image_indexmap .== 1] .= 0


cim = axs[1,1].imshow(df_f, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[1,1].set_title("ΔF/F")

cim = axs[2,1].imshow(med_filt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[2,1].set_title("""Remove "salt" noise""")

cim = axs[3,1].imshow(opened, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[3,1].set_title("Median filter")

# cim = axs[4,1].imshow(salt, cmap=cmap,
#     clim=(cmin,cmax), interpolation=nothing)
# axs[4,1].set_title("segments < 10")

cim = axs[4,1].imshow(df_f_nosalt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[4,1].set_title("Remove ROIs smaller than neurons")

for ax in axs
    ax.set_axis_off()
    outlines_to_show = small_outlines[:,:,z]
    outlines_to_show = Float64.(outlines_to_show)
    ax.imshow(outlines_to_show, alpha=outlines_to_show,cmap="binary_r")
end

scalebar = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="lower right", color = "white",
    pad=0.1, font_properties=Dict("size" => 7),
    scale_formatter = py"""lambda value, unit: "" """)
axs[1,1].add_artist(scalebar)

axs[1,1].text(3,40,"LSTM",fontsize=7, color="white")


plt.tight_layout()


analysis_name = "denoising_breakdown"
plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_df_f")
fig.savefig("$(plotpath).svg",
    dpi=1200)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=1200)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=1200)
@show plotpath*".svg"

fig
##

imshow(trial_average_raw[:,:,4,16,:])
imshow(trial_average_lstm[:,:,4,16,:])
size(trial_average_raw)

## temporal denoising
figB = Float64(uconvert(u"inch", 183u"mm") / u"inch")
figW = figB
figH = figB / 5 * 3
fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(3,5, figsize=(figW,figH))

timesteps = collect(1:8:35) .+ 2

for (j,t) in enumerate(timesteps)
    for (i,m) in enumerate([trial_average_raw[:,:,4,16,t], trial_average_kalman[:,:,4,16,t], trial_average_lstm[:,:,4,16,t]])
        cim = axs[i,j].imshow(m, interpolation=nothing, cmap="inferno")
        # axs[1,1].set_title("Raw")
        axs[i,j].set_axis_off()
    end
end
plt.tight_layout()

analysis_name = "denoising_temporal"
plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_df_f")

fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=300)
@show plotpath
fig
## temporal denoising - compare f0,f,df,df/f per raw, kalman, lstm
# FINAL FIG
window = 9
# window = window_len
figB = Float64(uconvert(u"inch", 183u"mm") / u"inch")
figW = figB
figH = figB / 4 * 3.2
fig = plt.figure(figsize=(figW,figH))
# nmodels x [f0, f, f-f0, df/f]
fig, axs = plt.subplots(3,4, figsize=(figW,figH))
ϵ=0.1
z = 3
for (i,(m,a,model_text)) in enumerate(zip(
    # [projective_field_raw[:,:,4,16], projective_field_kalman[:,:,4,16], projective_field_feng[:,:,4,16]],
    # [trial_average_raw[:,:,4,16,:], trial_average_kalman[:,:,4,16,:], trial_average_feng[:,:,4,16,:]])
    [projective_field_raw[:,:,z,16], projective_field_kalman[:,:,z,16], projective_field_lstm[:,:,z,16]],
    [trial_average_raw[:,:,z,16,:], trial_average_kalman[:,:,z,16,:], trial_average_lstm[:,:,z,16,:]],
    ["Raw", "Kalman", "LSTM"])
)
    f = mean(a[:,:,end-window+1:end],dims=3)[:,:,1]
    f0 = mean(a[:,:,1:window],dims=3)[:,:,1]
    df = f - f0
    # if i < 3
    if i < 3
        ϵ = 75
        ϵ = 0.0
        cmax4 = 2
    else
        ϵ = 0.0
        # cmax4 = maximum(m)
        cmax4 = 1
    end

    
    df_f = df./(f0 .+ ϵ)
    df_f_denoised = L.unsalt_n_pepper(df_f)
    cmax1 = maximum(f0)
    # cmax1 = quantile(f0[:],0.99)
    cmax2 = maximum(f)
    # cmax2 = quantile(f[:],0.99)
    # cmax3 = maximum(df)
    cmax3 = quantile(df[:],0.99)
    @show cmax1, cmax2, cmax3, cmax4, ϵ
    cim = axs[i,1].imshow(f0, clim=(0,cmax1), interpolation=nothing, cmap="inferno")
    axs[i,1].set_axis_off()
    cim = axs[i,2].imshow(f, clim=(0,cmax2), interpolation=nothing, cmap="inferno")
    axs[i,2].set_axis_off()
    cim = axs[i,3].imshow(df, clim=(0,cmax3), interpolation=nothing, cmap="inferno")
    axs[i,3].set_axis_off()
    cim = axs[i,4].imshow(df_f, clim=(0,cmax4), interpolation=nothing, cmap="inferno")
    axs[i,4].set_axis_off()
    # cim = axs[i,5].imshow(df_f_denoised, clim=(0,cmax4), interpolation=nothing, cmap="inferno")
    # axs[i,5].set_axis_off()
    for ax in [axs[i,j] for j in 1:4]
        outlines_to_show = small_outlines[:,:,z]
        outlines_to_show = Float64.(outlines_to_show)
        ax.imshow(outlines_to_show, alpha=outlines_to_show,cmap="binary_r")
    end

    axs[i,1].text(3,40,model_text,fontsize=7, color="white")
    if i == 1
        axs[1,1].set_title("F₀")
        axs[1,2].set_title("F")
        axs[1,3].set_title("ΔF = F - F₀")
        axs[1,4].set_title("ΔF/F")
        # axs[1,5].set_title("spatial denoising of ΔF/F")
        scalebar = matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
            scale_loc="left", location="lower right", color = "white",
            pad=0.1, font_properties=Dict("size" => 7),
            scale_formatter = py"""lambda value, unit: "" """)
        axs[1,1].add_artist(scalebar)
    end
end
plt.tight_layout()

analysis_name = "denoising_temporal_df_breakdown"
plotpath = joinpath(plot_dir,"$(analysis_name)_titration_oneplane_df_f")

fig.savefig("$(plotpath).svg",
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).png"),
    dpi=300)
fig.savefig(joinpath(plot_dir,"$(plotpath).pdf"),
    dpi=300)
@show plotpath
fig