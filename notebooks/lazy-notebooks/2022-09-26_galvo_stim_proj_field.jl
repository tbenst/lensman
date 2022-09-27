##
# ENV["DISPLAY"] = "localhost:11.0"
# ENV["DISPLAY"] = "/private/tmp/com.apple.launchd.5OQi0gJ6DL/org.xquartz:0"
##
# use for missing voltage / trigger data for using alternate galvo stim path
using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation,
    Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
    Unitful, StatsBase
using ImageView
import PyPlot
import Plots
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
inset_axes = pyimport("mpl_toolkits.axes_grid1.inset_locator").inset_axes
matscale = pyimport("matplotlib_scalebar.scalebar")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman

matplotlib.rcParams["font.sans-serif"] = ["Arial", "sans-serif"]

init_workers(16) # reduce memory usage
# init_workers(8) # reduce memory usage
resources = Resources();
##
r = Recordings[
    "2022-09-14_rschrmine_h2b6s_6dpf/fish3/TSeries_stim-thal_10z_1-aperature_4MHz-070" # no SLM trigger since using laser stim
](;
    resources...,
    tseries_read_strategy=:tiff
);

##
@pun (tseries,) = r;
##

# snippet to recover from missing or bad VoltageRecording
green_bot = mean(tseries, dims=(1, 2, 3))[1, 1, 1, :];
# red_bot = loadBOT(tseries_dir, "Ch2");
##
Plots.plot(green_bot)
##
thresh = 5 * std(green_bot)
putative_stim_starts = findall(green_bot .> thresh)
# not neighboring frames
valid_idx = vcat([1], diff(putative_stim_starts) .> 10) .== 1
putative_stim_starts = putative_stim_starts[valid_idx]
n_frames_between_stim = mode(diff(putative_stim_starts))
putative_stim_starts # no stim artifact for our control.... oops

##
# impute missing stimuli...
stim_start_idx = [putative_stim_starts[1]]
for i in 2:length(putative_stim_starts)
    isi = putative_stim_starts[i] - putative_stim_starts[i-1]
    prev = putative_stim_starts[i-1]
    # fill in missing stimuli
    while isi > 1.8 * n_frames_between_stim
        next = prev + n_frames_between_stim
        push!(stim_start_idx, next)
        isi -= n_frames_between_stim
        prev = next
    end
    push!(stim_start_idx, putative_stim_starts[i])
end


if length(trial_order) != length(stim_start_idx)
    # missing first control stims & last control stims...

    # beginning control stims...
    num_to_add = findfirst(trial_order .!= 3) - 1
    beg_start_idxs = []
    start_idx = stim_start_idx[1] - n_frames_between_stim
    while num_to_add > 0
        @info "beg" num_to_add, start_idx
        if start_idx > 0.7 * n_frames_between_stim
            pushfirst!(beg_start_idxs, start_idx)
        else
            @warn "start inference fails sanity check, ignoring" num_to_add, start_idx
        end
        start_idx -= n_frames_between_stim
        num_to_add -= 1
    end

    num_to_add = length(trial_order) - findlast(trial_order .!= 3)
    end_start_idxs = []
    start_idx = stim_start_idx[end] + n_frames_between_stim
    while num_to_add > 0
        push!(end_start_idxs, start_idx)
        start_idx += n_frames_between_stim
        num_to_add -= 1
    end
    stim_start_idx = vcat(beg_start_idxs, stim_start_idx, end_start_idxs)
end
stim_start_idx

@assert abs(maximum(diff(stim_start_idx)) - minimum(diff(stim_start_idx))) < 3