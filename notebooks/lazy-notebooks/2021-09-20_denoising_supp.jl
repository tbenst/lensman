ENV["DISPLAY"] = "localhost:10"
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
    trial_average_path="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_lstm_avgStim.h5",
    trial_average_dset="/lstm");
@pun (plot_dir, trial_average_path, trial_average_dset, tseries_units) = r_lstm
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
tseries_feng = r_feng[:tseries];
trial_average_lstm = r_lstm[:trial_average];
trial_average_raw = r_raw[:trial_average];
trial_average_kalman = r_kalman[:trial_average];
trial_average_feng = r_feng[:trial_average] .* 8192;
trial_average_feng = r_feng[:trial_average];
@show minimum(trial_average_feng), mean(trial_average_feng), maximum(trial_average_feng)
@show minimum(trial_average_raw), mean(trial_average_raw), maximum(trial_average_raw)
##
Gray.(imadjustintensity(tseries_feng[:,:,4,100]))
tseries_feng
##
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


# cmax = 0.02
# cmin = -0.02


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

## raw only
figB = Float64(uconvert(u"inch", 183u"mm") / u"inch")
figW = figB
figH = figB / 3
fig = plt.figure(figsize=(figW,figH))
fig, axs = plt.subplots(1,5, figsize=(figW,figH))
z = 4
s = 16
cmax = 1.0
cmin = 0
# cmin = -0.5
cmap = "viridis"
df_f = projective_field_raw[:,:,z,16]
# df_f = projective_field_kalman[:,:,z,16]
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
axs[1,1].set_title("Raw")

cim = axs[2,1].imshow(med_filt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[2,1].set_title("Morphological opening")

cim = axs[3,1].imshow(opened, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[3,1].set_title("Median filtered")

cim = axs[4,1].imshow(salt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[4,1].set_title("segments < 10")

cim = axs[5,1].imshow(df_f_nosalt, cmap=cmap,
    clim=(cmin,cmax), interpolation=nothing)
axs[5,1].set_title("spatial denoising")

for ax in axs
    ax.set_axis_off()
end
plt.tight_layout()


um_per_px1 = tseries_units[2] / 1u"μm"
scalebars = [matscale.ScaleBar(um_per_px1, "um", length_fraction=0.18, box_alpha=0,
    scale_loc="left", location="lower right", color = "white",
    font_properties=Dict("size" => 7),
    scale_formatter = py"""lambda value, unit: "" """)
    # )
for _ in 1:2]
axs[1].add_artist(scalebars[1])


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
