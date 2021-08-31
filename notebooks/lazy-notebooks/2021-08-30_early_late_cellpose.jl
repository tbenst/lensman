ENV["DISPLAY"] = "localhost:10"
##
using Lensman, PyCall, DataFrames, Gadfly, Distributed, StatsBase, Folds, ProgressMeter,
    Images, Statistics
import PyPlot
import Lensman: @pun, @assign
plt = PyPlot
np = pyimport("numpy")
L = Lensman
init_workers()
##
resources = Resources()
recording = Recordings["2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"](
    ; resources...
)
##

@pun (tseries_dir,
      trial_order, stim_start_idx, stim_end_idx, tseriesH, tseriesW,
      tseriesZ, vol_rate, nStimuli, nTrials, zseries, imaging2zseries_plane,
      zseriesH, zseriesW, tseries, tyh5_path, df_f_per_trial_dataframe,
      fish_dir, tseriesT
      ) = recording;
##
# cell_masks = L.read_cellpose_dir(joinpath(fish_dir, "cellpose_masks"));
cell_masks = L.read_cellpose_dir(joinpath(fish_dir, "cellpose_masks_conservative"));
cell_masks = shift_every_other_row(cell_masks,3); # undo from max z-proj
cell_ids = unique(cell_masks)
nCells = Int(length(cell_ids))
mask_idxs = Folds.map(c->findall(cell_masks .== c), cell_ids)
nCells = length(mask_idxs)
##
"Return nCell x T traces."
function cell_masks_to_traces(tseries::LazyTy5, mask_idxs)
    tyh5_path, dset_str = tseries.tyh5_path, tseries.dset_str
    nCells = length(mask_idxs)
    tseriesT = size(tseries,4)
    @showprogress @distributed hcat for t in 1:tseriesT
        proc_tseries = LazyTy5(tyh5_path, dset_str)
        vol = proc_tseries[:,:,:,t]
        neurons = [L.mean(vol[m]) for m in mask_idxs]
    end
end
neurons = cell_masks_to_traces(tseries, mask_idxs);
##

df_f_per_cell_per_trial = zeros(nTrials, nCells)
for (i,(st,en)) in enumerate(zip(stim_start_idx,stim_end_idx))
    window_len = Int(floor(5 * vol_rate)) - 1
    s = st - 1
    e = en + 1
    f0 = mean(neurons[:,s-window_len:s],dims=2)
    f = mean(neurons[:,e:e+window_len],dims=2)
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
threshs = collect(0:0.0001:0.005)
# tp = avg_stim_df[avg_stim_df[!,:stim] .== 16, :df_f_mean]
# fp = avg_stim_df[avg_stim_df[!,:stim] .== 1, :df_f_mean]
# threshs = collect(0.5:0.05:2.5)
fdrs = map(t->get_fdr(tp, fp, t), threshs)
# plot(x=threshs, y=fdrs, Geom.LineGeometry)
# nabove = map(t->sum(tp .> t)/length(tp), threshs)
nabove = map(t->sum(tp .> t), threshs)
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
    return -1
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
idxs = map(c-> c in recruited_ids, cell_ids)
recruited_cell_masks = mask_idxs[idxs];
# plot(recruited_stimthresh_df, x=:x1, Geom.histogram)

function label_by_cellid(mask_idxs, cell_ids,
        H, W, Z, value_lookup=identity)
    labels = ones(UInt32, H, W, Z)*-1
    
    for (c,mask) in zip(cell_ids, mask_idxs)
        labels[mask] .= value_lookup(c)
    end
    labels
end

recruited_labels = label_by_cellid(recruited_cell_masks, recruited_ids,
    tseriesH, tseriesW, tseriesZ,
    c->getindex(Dict(eachrow(recruited_stimthresh_df[:,[:cell_id, :x1]])) ,c));

     
##
tseries_subset = tseries[:,:,:,1:20:tseriesT];
maxTproj = mapslices(x->quantile(x[:],0.90), tseries_subset,dims=4);
shiftedmaxTproj = shift_every_other_row(maxTproj,-3);
shiftedmaxTproj = shiftedmaxTproj[:,:,:,1];


##
fig, ax = plt.subplots(1,dpi=600)
p = 7
im = maximum(shiftedmaxTproj[175:310,100:235,:], dims=3)[:,:,1]
im = L.adjust_histogram(im, GammaCorrection(0.2))
# ax.imshow(zseries[:,:,imaging2zseries_plane[p]], cmap=plt.cm.gray, alpha=0.2)
ax.imshow(im, cmap=plt.cm.gray)
# ax.imshow(maximum(zseries[:,:,imaging2zseries_plane], dims=3),
#     cmap=plt.cm.gray)
# ax.imshow(imresize(recruited_labels[:,:,p],zseriesH, zseriesW),
#     cmap=plt.cm.magma_r, alpha=0.3)
# rl = recruited_labels
rl = shift_every_other_row(recruited_labels,-3)
rl = maximum(rl, dims=3)[175:310,100:235,1] .* 2
rl = pycall(np.ma.masked_where, Any, rl .== -2, rl)

cim = ax.imshow(rl, cmap=plt.cm.magma_r, alpha=0.6)
cbar = plt.colorbar(cim)
cbar.set_label("# of stimulated targets")
ax.set_xticks([])
ax.set_xticks([], minor=true)
ax.set_yticks([])
ax.set_yticks([], minor=true)
fig.savefig("/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/2021-06-08_rsChRmine_titration/cell-stim-threshold_cellposed.png")
fig

##
x = zeros(size(tseries)[1:3])
for m in mask_idxs[2:end]
    x[m] .= 1
end
Gray.(maximum(x,dims=3)[:,:,1])