# need X11 forwarding for ssh, `ssh -X osprey`, then `echo $DISPLAY`, and fill below
ENV["DISPLAY"] = "localhost:10.0"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, HDF5,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman
using Plots
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
matplotlib.use("tkagg")
using Base.Iterators:peel
import Unitful:μm

## Define experiment paths

# November 3-region stim
fish_dir = "/scratch/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2"
exp_name = "TSeries-lrhab_raphe_40trial-045"

# More recent 3-region stim with 4th control region
# fish_dir = "/scratch/b115/2021-02-09_gcamp6f_7dpf/fish2"
# exp_name = "TSeries-lrhab-raphe-balanced-transitions-137"

tif_dir = joinpath(fish_dir, exp_name)

## Load tseries
tseries = loadTseries(tif_dir)
(H, W, Z, T) = size(tseries)

## Get stim times from voltage recordings
voltage_file = glob("*VoltageRecording*.csv", tif_dir)[1]
stim_start_idx, stim_end_idx = getStimTimesFromVoltages(voltage_file, Z)

## Parse experiment metadata
data_folders = splitpath(tif_dir)

# e.g. <fish_dir>/<exp_name>/<exp_name>.xml 
xml_path = joinpath(data_folders..., data_folders[end] * ".xml")
exp_date, frame_rate, etl_vals = getExpData(xml_path)

# Volume gets imaged at 1/Zth of the frame rate
vol_rate = frame_rate / Z

# if imaging many planes, may need to read another xml file since didn't get all planes in first file
@assert length(etl_vals) == Z
@info "assume imaging from top-down"

## Load SLM files
# slm_dir = "/mnt/deissero/users/tyler/b115/SLM_files"
slm_dir = "/scratch/allan/b115"
slm_exp_dir = joinpath(slm_dir, Dates.format(exp_date, "dd-u-Y"))

trial_order, _ = getTrialOrder(slm_exp_dir, exp_date)  # in one case reassigns slm_exp_dir

n_stimuli = maximum(trial_order)
n_trials = size(trial_order, 1)
n_trials_per_stimulus = Int(size(trial_order, 1) / n_stimuli)
@assert n_trials == size(stim_start_idx, 1)







## Explore trial_order to see how many of each sequence we have

stim_type_to_name = Dict([(1, "left_habenula"), (2, "right_habenula"), (3, "raphe")])

two_stim_sequences = [(first, second) for first = 1:n_stimuli for second = 1:n_stimuli]

# Lookup from sequence of two stim types (e.g. (1, 1)) to timestep in which second stim started
# (implicitly we are discarding first stim in trial as there was no other stim that came before it)
sequence_to_stim_index = Dict([(x, Int[]) for x in two_stim_sequences])
for stim_number in 1:length(trial_order) - 1
    key = (trial_order[stim_number], trial_order[stim_number + 1])
    push!(sequence_to_stim_index[key], stim_start_idx[stim_number + 1])
end


## High-level viz settings
cmax = 8
cmin = -0.9
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)

# Size of window over which we are computing mean f and f0
pre = Int(ceil(5 * vol_rate))
post = pre
window = Int(ceil(3 * vol_rate))
@assert (window < post) & (window < pre)



##
# Standard baseline fluorescence from before
# avg_stim_per_trial_type = trialAverage(tseries, stim_start_idx, stim_end_idx, trial_order;
#     pre=pre, post=post);



## Trial average for each sequence
# TODO: remove
pre = 16
post = 50
sequence_to_trial_average = Dict()
for sequence in two_stim_sequences
    sequence_to_trial_average[sequence] = trial_average_select(
        tseries, stim_start_idx, stim_end_idx, sequence_to_stim_index[sequence], pre=pre, post=post)
end
for (k, v) in sequence_to_trial_average
    println("$k: $(size(v))")
end

## Make a video
# plane = 1
# imgstack = []
# sequence = (3, 2) # (1, 2)
# trial_average = sequence_to_trial_average[sequence][1:end, 1:end, plane, 1:end]
# trial_average = floor.(UInt8, trial_average / maximum(trial_average) * 255)
# println(typeof(trial_average))
# println("$(minimum(trial_average)) $(maximum(trial_average))")
# using VideoIO
# for _t in 1:size(trial_average)[3]
#     image = trial_average[1:end, 1:end, _t]
#     println(typeof(image))
#     push!(imgstack, image)
# end
# props = [:color_range => 2, :priv_data => ("crf" => "0", "preset" => "medium")]
# println(length(imgstack))
# encodevideo("$(stim_type_to_name[sequence[2]])_after_$(stim_type_to_name[sequence[1]]).mp4", imgstack, framerate=3, AVCodecContextProperties=props)

plane = 1
imgstack = []
sequences = [(1, 2), (2, 2), (3, 2)] # (1, 2)
# trial_average = sequence_to_trial_average[sequence][1:end, 1:end, plane, 1:end]
# trial_average = floor.(UInt8, trial_average / maximum(trial_average) * 255)
# println(typeof(trial_average))
# println("$(minimum(trial_average)) $(maximum(trial_average))")
# maxi = maximum([maximum(sequence_to_trial_average[_s]) for _s in sequences])
maxi = 1000
print(maxi)
using VideoIO
nt = size(sequence_to_trial_average[sequences[1]])[4]
for _t in 1:nt
    _seqs = [sequence_to_trial_average[_s][1:end, 1:end, plane, _t] for _s in sequences]
    image = hcat(_seqs[1], _seqs[2], _seqs[3])
    println(size(image))
    image = floor.(UInt8, clamp.(image, 0, 1000) / maxi * 255)
    println(typeof(image))
    push!(imgstack, image)
end
props = [:color_range => 2, :priv_data => ("crf" => "0", "preset" => "medium")]
println(length(imgstack))
# encodevideo("$(stim_type_to_name[sequence[2]])_after_$(stim_type_to_name[sequence[1]]).mp4", imgstack, framerate=3, AVCodecContextProperties=props)
encodevideo("right_habenula.mp4", imgstack, framerate=3, AVCodecContextProperties=props)

## First just coarse viz (with one influence map for each stim sequence)
cmax = 8
cmin = -0.9
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)
# fig = plt.figure(figsize=(50, 30))
fig = plt.figure()
for (idx, sequence) in enumerate(two_stim_sequences)

    avg_stim = sequence_to_trial_average[sequence]

    # Looking at a specific window post stim (NOTE: this really would be more apples to apples,
    # but likely a difference of just a couple of frames)
    # f = mean(avg_stim[:, :, :, end - window + 1:end], dims=4)[:, :, :, 1]

    # Instead want to look at first index post stim onwards
    f = mean(avg_stim[:, :, :, pre + 8:end], dims=4)[:, :, :, 1]   # approx

    # Baseline computed as all trials of this type WITH previous type as first type in sequence
    f0 = mean(avg_stim[:, :, :, 1:window], dims=4)[:, :, :,1]

    # Baseline fluorescence computed across all trials of this type
    # f0 = mean(avg_stim_per_trial_type[:, :, :, sequence[2], 1:window], dims=4)[:, :, :,1]

    df = f - f0
    df_f = df ./ f0

    cmax = percentile(df_f[:], 99.9)
    cmin = percentile(df_f[:], 0.1)
    println(cmax, cmin)

    # plt.subplot(length(two_stim_sequences), 1, idx)
    plt.axis("off")
    plt.imshow(hcat([df_f[:, :, z] for z in 1:Z]...), cmap="RdBu_r", norm=cnorm)
    plt.title("$(stim_type_to_name[sequence[2]]) stim (after $(stim_type_to_name[sequence[1]]) stim)")
    fig.savefig("viz/$(stim_type_to_name[sequence[2]])_after_$(stim_type_to_name[sequence[1]])", bbox_inches="tight", pad_inches=0)
end
# display(fig)





## Trial average (NOTE: these are dropping some frames right after stim)

pre = 16
post = 16
avg_stim = trialAverage(tseries, stim_start_idx, stim_end_idx, trial_order;
    pre=pre, post=post);

## Coarse effect
window = Int(ceil(3 * vol_rate))  # 10
cmax = 4
cmin = -0.8
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)

fig = plt.figure(figsize=(50, 30))
for stim_index in 1:n_stimuli

    f0 = mean(avg_stim[:, :, :, stim_index, 1:window], dims=4)[:, :, :,1]
    f = mean(avg_stim[:, :, :, stim_index, end - window + 1:end], dims=4)[:, :, :, 1]
    df = f - f0
    df_f = df ./ f0

    cmax = percentile(df_f[:], 99.9)
    cmin = percentile(df_f[:], 0.1)
    println(cmax, cmin)

    plt.subplot(3, 1, stim_index)
    plt.axis("off")
    plt.imshow(hcat([df_f[:, :, z] for z in 1:Z]...), cmap="RdBu_r",
        norm=cnorm)

end
display(fig)
fig.savefig("stim_coarse_average.png", bbox_inches="tight", pad_inches=0)


##

# High-level viz settings
cmax = 8
cmin = -0.9
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)

# Size of window over which we are computing mean f and f0
window = Int(ceil(3 * vol_rate))
@assert (window < post) & (window < pre)

# Length of block from 5s before stim, through stim, until 5s after stim (if post=16
# otherwise can make up to next stim)
block_length = size(avg_stim)[end]

# Duration of a stim in timesteps
stim_length = minimum(stim_end_idx - stim_start_idx) + 1

# Parameters
# first_index = pre + stim_length + 5  # avoids stim artifact
# n_indices = block_length - first_index + 1
# n_blocks = 10

# Trying now:
first_index = pre + stim_length + 1
n_indices = block_length - first_index + 1
n_blocks = n_indices - 1

block_size = Int(floor(n_indices / n_blocks))

indices = collect(first_index:block_size:block_length)
println("Block size: $block_size, num_bins: $(length(indices))")
println(first_index)
println(indices)
##

fig = plt.figure(figsize=(50, 30))
for stim_index in 1:n_stimuli

    f0 = mean(avg_stim[:, :, :, stim_index, 1:window], dims=4)[:, :, :,1]
    
    for bin in 1:n_blocks - 1
        start_index = indices[bin]
        end_index = indices[bin + 1]

        f = mean(avg_stim[:, :, :, stim_index, start_index:end_index], dims=4)[:, :, :, 1]

        df = f - f0
        df_f = df ./ f0

        cmax = percentile(df_f[:], 99.9)
        cmin = percentile(df_f[:], 0.1)
        println(cmax, cmin)

        # plt.subplot(n_blocks, 1, stim_index)
        plt.subplot(1, n_blocks, bin)
        plt.axis("off")
        plt.imshow(vcat([df_f[:, :, z] for z in 1:Z]...), cmap="RdBu_r",
            norm=cnorm)

    end

    # fig.savefig("stim_$stim_index.png")
    display(fig)
end






## Load target groups
target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slm_exp_dir))]
is1024 = true
targets_with_plane_index = mapTargetGroupsToPlane(target_groups, etl_vals,
    is1024=is1024, zOffset=0)

## Visualize cells
@warn "Hard-coding microscope units for 512x512"
# TODO: ahh important to fix
microscope_units_xy = 0.6299544139175637μm
if is1024
    microscope_units_xy *= 2
end
# assume that "7um" spiral galvo signal is calibrated for 16x, not 25x
# "16x" is actually 14.4x
target_size_px = 7μm * (14.4 / 25) / microscope_units_xy

# Generate mask of targeted locations
stim_masks = Gray.(zeros(Bool, H, W, Z, n_stimuli))
for (stim_idx, target_group) in enumerate(targets_with_plane_index)
    for (x, y, z) in eachrow(target_group)
        draw!(view(stim_masks, :, :, z, stim_idx), Ellipse(CirclePointRadius(x, y, target_size_px)))
    end
end

##
window = Int(ceil(3 * vol_rate))
pre = 16
post = minimum(stim_start_idx[2:end] - stim_start_idx[1:end - 1]) - 10  # comes out to 50

# Generate masks
stim_masks_one_plane = zeros(Bool, H, W, Z, n_stimuli)
stim_masks_one_plane[:, :, :, 1] = stim_masks[:, :, :, 1]
stim_masks_one_plane[:, :, :, 2] = stim_masks[:, :, :, 2]
stim_masks_one_plane[:, :, :, 3] = stim_masks[:, :, :, 3]


stim_start_to_end_lookup = Dict(stim_start_idx .=> stim_end_idx)

for (idx, sequence) in enumerate(two_stim_sequences)
    indices = sequence_to_stim_index[sequence]

    region_to_F_sequence = Dict([(region, []) for region in 1:3])

    for start in indices
        first_index_post_stim = start + (stim_start_to_end_lookup[start] - start) + 1

        for region in 1:3


            # # First approach, not bad but noisy:
            # f0 = mean(stim_masks_one_plane[:, :, :, region] .* view(tseries, :, :, :, start - window:start - 1))
            # influence = stim_masks_one_plane[:, :, :, region] .* view(tseries, :, :, :, first_index_post_stim:first_index_post_stim + post - 1)
            # mn = (mean(influence, dims=[1,2,3])[1, 1, 1, :] .- f0) ./ f0
            # push!(region_to_F_sequence[region], mn)


            # # First produce 1D time series then take Kalman filter of that
            # look_back = 5
            # tseries_1d = view(tseries, :, :, :, start - window - look_back:start - 1) .* stim_masks_one_plane[:, :, :, region]
            # tseries_1d = mean(tseries_1d, dims=[1,2,3])[1, 1, 1, :]
            # baseline_kf = imageJkalmanFilter(tseries_1d)[look_back + 1:end]
            # f0 = mean(baseline_kf)
            # # Running KF on subset of tseries
            # tseries_1d = view(tseries, :, :, :, first_index_post_stim - look_back:first_index_post_stim + post - 1) .* stim_masks_one_plane[:, :, :, region]
            # tseries_1d = mean(tseries_1d, dims=[1,2,3])[1, 1, 1, :]
            # tseries_kf = imageJkalmanFilter(tseries_1d)[look_back + 1:end]
            # mn = (tseries_kf .- f0) ./ f0
            # push!(region_to_F_sequence[region], mn)

            

            # look_back = 5  # if this is too big will be looking at previous stim which is bad
            # num_ones_in_mask = sum(stim_masks_one_plane[:, :, :, region])
            # tseries_1d = view(tseries, :, :, :, start - window - look_back:start - 1) .* stim_masks_one_plane[:, :, :, region]
            # tseries_1d = sum(tseries_1d, dims=[1,2,3])[1, 1, 1, :] ./ num_ones_in_mask
            # baseline_kf = imageJkalmanFilter(tseries_1d)[look_back + 1:end]


            # # Trying to just plot Kalman filtered fluorescence  (NOTE: going into stim region with burn in)
            # # tseries_1d = view(tseries, :, :, :, first_index_post_stim - look_back:first_index_post_stim + post - 1) .* stim_masks_one_plane[:, :, :, region]
            # tseries_1d = view(tseries, :, :, :, first_index_post_stim:first_index_post_stim + post - 1) .* stim_masks_one_plane[:, :, :, region]
            # tseries_1d = sum(tseries_1d, dims=[1,2,3])[1, 1, 1, :] ./ num_ones_in_mask  # Ok to mean and then filter?
            # # tseries_kf = imageJkalmanFilter(tseries_1d)[look_back + 1:end]
            # tseries_kf = imageJkalmanFilter(tseries_1d)

            # # push!(region_to_F_sequence[region], tseries_kf)
            # push!(region_to_F_sequence[region], vcat(baseline_kf, tseries_kf))


            # This could introduce stim artifacts
            plot_start = start - pre
            plot_end = stim_start_to_end_lookup[start] + post
            plot_range = plot_start:plot_end
            trace = view(tseries, :, :, :, plot_range) .* stim_masks_one_plane[:, :, :, region]
            num_ones_in_mask = sum(stim_masks_one_plane[:, :, :, region])
            trace = sum(trace, dims=[1,2,3])[1, 1, 1, :] ./ num_ones_in_mask
            # trace = extractTrace(tseries[:, :, :, plot_range], stim_masks_one_plane[:, :, :, region])
            trace = imageJkalmanFilter(trace)
            print(length(trace))
            push!(region_to_F_sequence[region], trace[1:73])  # hacky but should work
        end
    end
    # p = plot(plot_range / vol_rate, trace, left_margin=50px, ylims=(0, 500))

    fig = plt.figure(figsize=(50, 10))
    for region in 1:3
        plt.subplot(1, 3, region)
        plt.ylim(50, 300)
        for mns in region_to_F_sequence[region]
            plt.plot((1:73) / vol_rate, mns)
        end
        plt.plot((1:73) / vol_rate, mean(region_to_F_sequence[region]), color=:black)

        # TODO: can always add rectangle later
        # plot!(Lensman.rectangle(50, 250, 16, 20), opacity=.5, label="stim")

        plt.title("$(stim_type_to_name[sequence[2]]) stim (after $(stim_type_to_name[sequence[1]]) stim), F in $(stim_type_to_name[region])")
    end
    fig.savefig("viz/$(stim_type_to_name[sequence[2]])_after_$(stim_type_to_name[sequence[1]])_trace", bbox_inches="tight", pad_inches=0)

end






##
# Hav


# Need to make sure that post never takes us into the next stimulation
# @assert n_trials == length(stim_start_idx)
# for stim_index in 2:n_trials - 1
#     stim_start = stim_start_idx[stim_index]
#     stim_end = stim_end_idx[stim_index]

#     # Using maximum to be conservative here (`trialAverage` uses minimum)
#     stim_length = maximum(stim_end_idx - stim_start_idx)

#     # Make sure we're not including next stim in our average (or previous stim in baseline)
#     # (second assert keeps a distance from end of previous stim to be careful about including
#     # remaining decay from previous stim)
#     @assert stim_start + stim_length + post < stim_start_idx[stim_index + 1]
#     @assert stim_start - pre - stim_end_idx[stim_index - 1] >= 37
# end



            # This produces issues with 0's
            # f0 = mean(stim_masks_one_plane[:, :, :, region] .* view(tseries, :, :, :, start - window:start - 1), dims=4)[:, :, :, 1]
            # println("Size of mask: $(sum(stim_masks_one_plane[:, :, :, region]))")
            # println("numzeros: $(sum(f0 .== 0) - sum(1 .- stim_masks_one_plane[:, :, :, region]))")
            # influence = stim_masks_one_plane[:, :, :, region] .* view(tseries, :, :, :, first_index_post_stim:first_index_post_stim + post - 1)
            # df = influence .- f0
            # println("df: $(size(df))")
