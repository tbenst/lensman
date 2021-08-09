using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics
import Lensman: @pun, @assign
import Plots
import Plots: heatmap
using Thunks
np = pyimport("numpy")
L = Lensman
# using Dagger
# import Dagger: @par
# import Dagger.Sch: ThunkOptions
##

r = Recording(
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123";
    zseries_name="ZSeries-structural-840nm-058",
    # tseries_dset = "/imaging/LSTM_per-voxel-state_divide2048-2021-07-02",
    tseries_dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm",
    tyh5_path="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-1232021-06-21_6pm.ty.h5",
    suite2p_dir="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/s2p/suite2p"
    # tseries_root_dirs = [
    #     "/data/b115"
    # ],
);
@pun (imaging2zseries_plane,) = r
##
@pun (nCells, cell_masks, tseries_dir, suite2p_dir, cell_traces,
      trial_order, stimStartIdx, stimEndIdx, tseriesH, tseriesW,
      tseriesZ, vol_rate, nStimuli, nTrials, zseries, imaging2zseries_plane
      ) = r;

## pattern for extending DAG (ugly! ergonomics should be improved..macro?)
if ~in(:vol_rate, keys(r.thunks))
    @thunk begin
        vol_rate = frame_rate / tseriesZ
    end
    @assign r.thunks = vol_rate
end
@pun        vol_rate = r

##
avg_tseries = avgTseries(tseries_dir);


##
neuron_labels = zeros(UInt16, tseriesH, tseriesW, tseriesZ)
# TODO: find non unique XYZ tuples & keep row with largest weight.

for c in unique(cell_masks.cell_id)
    pix_mat = Array(cell_masks[cell_masks.cell_id .== c,[:y,:x,:z]]) .+ 1
    voxel_mask = mapslices(x->CartesianIndex(x...),pix_mat,dims=2)
    neuron_labels[voxel_mask] .= c
end
##
colored_neurons = map(get_random_color, (neuron_labels[:,:,6]));
# colored_neurons
shift_every_other_row(colored_neurons,-3)

##
# ntime x nCells
size(cell_traces)
df_f_per_cell_per_trial = zeros(nTrials, nCells)
for (i,(st,en)) in enumerate(zip(stimStartIdx,stimEndIdx))
    window_len = Int(floor(5 * vol_rate)) - 1
    s = st - 1
    e = en + 1
    f0 = mean(cell_traces[s-window_len:s,:],dims=1)
    f = mean(cell_traces[e:e+window_len,:],dims=1)
    df_f_per_cell_per_trial[i,:] = @. (f - f0) / (f0 + 5)
end
df_f_per_cell_per_trial

df = DataFrame(df_f_per_cell_per_trial,
            map(string,collect(1:nCells)))
df[!, :stim] = trial_order
df = stack(df, 1:nCells)
rename!(df, Dict(:variable => :cell_id, :value => :df_f))
df[!,:cell_id] = parse.(Int,df[!,:cell_id])
@show size(df)
first(df,20)
##
avg_stim_df = combine(groupby(df, [:cell_id, :stim]), :df_f => mean)
idxs = (avg_stim_df.stim .== 1) .| (avg_stim_df.stim .== 16)
plot(avg_stim_df[idxs,:], x=:df_f_mean, ygroup=:stim, Scale.y_log10,
    Geom.subplot_grid(Geom.histogram))
##
min_stim_df = combine(groupby(df, [:cell_id, :stim]), :df_f => minimum)
idxs = (min_stim_df.stim .== 1) .| (min_stim_df.stim .== 16)
plot(min_stim_df[idxs,:], x=:df_f_minimum, ygroup=:stim, Scale.y_log10,
    Geom.subplot_grid(Geom.histogram))
##
function get_fdr(true_positives::Vector{<:Real}, false_positives::Vector{<:Real},
        threshold::Real)
    ntrue_positives = sum(true_positives .> threshold)
    nfalse_positives = sum(false_positives .> threshold)
    fdr = nfalse_positives / (ntrue_positives + nfalse_positives)
end

#TODO maybe: 0.2 quantile?

tp = min_stim_df[min_stim_df[!,:stim] .== 16, :df_f_minimum]
fp = min_stim_df[min_stim_df[!,:stim] .== 1, :df_f_minimum]
threshs = collect(0:0.01:0.5)
# tp = avg_stim_df[avg_stim_df[!,:stim] .== 16, :df_f_mean]
# fp = avg_stim_df[avg_stim_df[!,:stim] .== 1, :df_f_mean]
# threshs = collect(0.5:0.05:2.5)
fdrs = map(t->get_fdr(tp, fp, t), threshs)
# plot(x=threshs, y=fdrs, Geom.LineGeometry)
nabove = map(t->sum(tp .> t)/length(tp), threshs)
p1 = plot(x=threshs, y=fdrs, Geom.LineGeometry, Guide.title("fdr"))
p2 = plot(x=threshs, y=nabove, Geom.LineGeometry, Guide.title("nabove"))
vstack(p1,p2)


##
"Given a threshould, find the lowest stim number where all df_f are above"
function stim_threshold(cell_subdf, key, thresh, nStimuli)
    cell_id = cell_subdf[1,:cell_id]
    for s in 1:nStimuli
        idxs = cell_subdf.stim .>= s
        # if all(cell_subdf[idxs,key] .> thresh)
        if (sum(cell_subdf[idxs,key] .> thresh)/length(cell_subdf[idxs,key])) > 0.75
            return s
        end
    end
    # hack so background is black with 16 stimulation conditions
    return 17
end

thresh = 0.00
stimthresh_df = combine(c->stim_threshold(c,:df_f_minimum, thresh, nStimuli),
    groupby(min_stim_df, :cell_id))
# thresh = 1.00
# stimthresh_df = combine(c->stim_threshold(c,:df_f_mean, thresh, nStimuli),
#     groupby(avg_stim_df, :cell_id))

sum((~).(isnothing.(stimthresh_df.x1)))
idxs = stimthresh_df.x1 .>= 0
recruited_stimthresh_df = stimthresh_df[idxs,:]
recruited_ids = unique(recruited_stimthresh_df.cell_id)
idxs = map(c-> c in recruited_ids, cell_masks.cell_id)
recruited_cell_masks = cell_masks[idxs,:];
# plot(recruited_stimthresh_df, x=:x1, Geom.histogram)

function label_by_cellid(cell_masks::DataFrame,
        H, W, Z, value_lookup=identity)
    labels = ones(UInt16, H, W, Z)*17
    
    for c in unique(cell_masks.cell_id)
        pix_mat = Array(cell_masks[cell_masks.cell_id .== c,[:y,:x,:z]]) .+ 1
        voxel_mask = mapslices(x->CartesianIndex(x...),pix_mat,dims=2)
        labels[voxel_mask] .= value_lookup(c)
    end
    labels
end

recruited_labels = label_by_cellid(recruited_cell_masks,tseriesH, tseriesW, tseriesZ,
    c->getindex(Dict(eachrow(recruited_stimthresh_df[:,[:cell_id, :x1]])) ,c));

# Geom.histogram2d

# recruited_neurons = map(get_random_color,
#     maximum(recruited_labels,dims=3)[:,:,1])

# shift_every_other_row(colored_neurons,-3)
##
# frame_df = DataFrame(recruited_labels[:,:,6],string.(collect(1:tseriesW)))
# frame_df.y = collect(1:tseriesH)
# frame_df = stack(frame_df,1:tseriesW)
# rename!(frame_df, :variable => :x)
# plot(frame_df, x=:x, y=:y, Geom.rectbin, Coord.cartesian(raster=true))
##

heatmap(minimum(shift_every_other_row(recruited_labels,-2),dims=3)[150:350,100:220,1],
    c=Plots.cgrad(:thermal, rev = true))
# heatmap(recruited_labels[:,:,5])
# threshold of 1-1.4 seems reasonable
##

## example for extending DAG....
# @thunk begin
#     suite2p_dir = L.is_setting(:suite2p_dir, get_suite2p_dir(r.thunks[:tseries_dir]))
# end
##
suite2p_dir = thunk(L.is_setting)(r.settings, :suite2p_dir, thunk(L.get_suite2p_dir)(r.thunks[:tseries_dir]))
@assign r.thunks = suite2p_dir
@pun suite2p_dir = r

