ENV["DISPLAY"] = "localhost:11"
##
using ImageView
using Lensman, PyCall, DataFrames, Gadfly, Distributed, StatsBase, Folds, ProgressMeter,
    Images, Statistics, CategoricalArrays
import PyPlot
import Lensman: @pun, @assign
import Base: >, <, ==
plt = PyPlot
matplotlib = plt.matplotlib
np = pyimport("numpy")
L = Lensman
init_workers()
##
resources = Resources()
recording = Recordings["2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"](
    ; resources...,
    # tseries_read_strategy = :lazy_hwzt,
    # tyh5_path="/data/dlab/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123_kalman.h5",
    # tseries_dset=nothing

)
##
# Let's tak ea single ROI where I can see by eye what's happening
# debug code for a single-cell
# question: can we replicate first fig with this code?
# T
##

@pun (tseries_dir,
      trial_order, stim_start_idx, stim_end_idx, tseriesH, tseriesW,
      tseriesZ, vol_rate, nStimuli, nTrials, zseries, imaging2zseries_plane,
      zseriesH, zseriesW, tseries, tyh5_path, fish_dir, tseriesT, window_len, 
      plot_dir
) = recording;
##
# cell_masks = L.read_cellpose_dir(joinpath(fish_dir, "cellpose_masks"));
# cell_masks = L.read_cellpose_dir(joinpath(fish_dir, "cellpose_masks_conservative"));
cell_masks = L.read_cellpose_dir(joinpath(fish_dir, "cellpose_masks_kalman_v3"));
cell_masks = shift_every_other_row(cell_masks,3); # undo from max z-proj
cell_ids = unique(cell_masks)
nCells = Int(length(cell_ids))
mask_idxs = Folds.map(c->findall(cell_masks .== c), cell_ids)
big_enough = length.(mask_idxs) .> 9
mask_idxs = mask_idxs[big_enough]
cell_ids = cell_ids[big_enough]
nCells = length(mask_idxs)
##
"Return nCell x T traces."
# function cell_masks_to_traces(tseries::LazyTy5, mask_idxs)
function cell_masks_to_traces(tseries::Lensman.LazyHDF5, mask_idxs)
    tyh5_path, dset_str = tseries.tyh5_path, tseries.dset_str
    nCells = length(mask_idxs)
    tseriesT = size(tseries,4)
    @showprogress @distributed hcat for t in 1:tseriesT
        proc_tseries = typeof(tseries)(tyh5_path, dset_str)
        vol = proc_tseries[:,:,:,t]
        neurons = [L.mean(vol[m]) for m in mask_idxs]
    end
end
neurons = cell_masks_to_traces(tseries, mask_idxs);
# neurons_f0 = mapslices(x->quantile(x,0.1), neurons; dims=2)[:,1];
##
neurons
@warn "actually df not df/f"
df_f_per_cell_per_trial = zeros(nTrials, nCells)
for (i,(st,en)) in enumerate(zip(stim_start_idx,stim_end_idx))
    window_len = Int(floor(2 * vol_rate)) - 1
    # s = st - Int(floor(5 * vol_rate))
    # e = en - Int(floor(5 * vol_rate))
    s = st - 1
    e = en + 1
    f0 = mean(neurons[:,s-window_len:s],dims=2)
    f = mean(neurons[:,e:e+window_len],dims=2)
    # df_f_per_cell_per_trial[i,:] = @. (f - f0) / (f0 + 1-8)
    df_f_per_cell_per_trial[i,:] = (f - f0)
end
df_f_per_cell_per_trial

df_f_per_cell_per_trial = DataFrame(df_f_per_cell_per_trial,
            map(string,collect(1:nCells)))
df_f_per_cell_per_trial[!, :stim] = trial_order
df_f_per_cell_per_trial = stack(df_f_per_cell_per_trial, 1:nCells)
rename!(df_f_per_cell_per_trial, Dict(:variable => :cell_id, :value => :df_f))
df_f_per_cell_per_trial[!,:cell_id] = parse.(Int,df_f_per_cell_per_trial[!,:cell_id])
@show size(df_f_per_cell_per_trial)
first(df_f_per_cell_per_trial,20)

function add_period_to_df(df)
    df[!,:period] = ["" for _ in 1:size(df,1)];
    @assert maximum(df[!,:stim]) == 16
    for s in 1:16
        if s <= 6
            p = "early"
        elseif s <= 11
            p = "mid"
        else
            p = "late"
        end
        idxs = df[!, :stim] .== s
        df[idxs,:period] .= p
    end
    catarr = CategoricalArray(df[!,:period], ordered=true)
    levels!(catarr, ["early", "mid", "late"])
    df[!,:period] = catarr
    df
end

df_f_per_cell_per_trial = add_period_to_df(df_f_per_cell_per_trial);

avg_stim_df = combine(groupby(df_f_per_cell_per_trial, [:cell_id, :stim]), :df_f => mean =>:df_f)
avg_stim_per_period_df = combine(groupby(df_f_per_cell_per_trial, [:cell_id, :period]), :df_f => mean =>:df_f)

avg_stim_df = add_period_to_df(avg_stim_df)
min_stim_df = combine(groupby(df_f_per_cell_per_trial, [:cell_id, :stim]), :df_f => minimum =>:df_f)
min_stim_df = add_period_to_df(min_stim_df)
med_stim_df = combine(groupby(df_f_per_cell_per_trial, [:cell_id, :stim]), :df_f => median =>:df_f)
med_stim_df = add_period_to_df(med_stim_df)
mean(df.df_f)
plot(med_stim_df, y=:df_f, x=:period, Geom.boxplot(suppress_outliers=true),
    Coord.Cartesian(;ymin=-0.001,ymax=0.001))
plot(avg_stim_df, y=:df_f, x=:period, Geom.boxplot(suppress_outliers=true),
    Coord.Cartesian(;ymin=-0.001,ymax=0.001))
plot(avg_stim_per_period_df, y=:df_f, x=:period, Geom.boxplot(suppress_outliers=true),
    Coord.Cartesian(;ymin=-0.001,ymax=0.001))
##
function get_fdr(true_positives::Vector{<:Real}, false_positives::Vector{<:Real},
        threshold::Real)
    ntrue_positives = sum(true_positives .> threshold)
    nfalse_positives = sum(false_positives .> threshold)
    fdr = nfalse_positives / (ntrue_positives + nfalse_positives)
end

#TODO maybe: 0.2 quantile?

threshs = collect(0:0.00001:0.0005)
the_stim_df = med_stim_df
tp = the_stim_df[the_stim_df[!,:stim] .== 16, :df_f]
fp = the_stim_df[the_stim_df[!,:stim] .== 1, :df_f]
# threshs = collect(0.5:0.05:2.5)
fdrs = map(t->get_fdr(tp, fp, t), threshs)
# plot(x=threshs, y=fdrs, Geom.LineGeometry)
# nabove = map(t->sum(tp .> t)/length(tp), threshs)
nabove = map(t->sum(tp .> t), threshs) / length(tp)
p1 = plot(x=threshs, y=fdrs, Geom.LineGeometry, Guide.title("fdr"))
p2 = plot(x=threshs, y=nabove, Geom.LineGeometry, Guide.title("nabove"))
vstack(p1,p2)

##
"Given a threshould, find the lowest stim number where all df_f are above"
function stim_threshold(cell_subdf, key, thresh; frac=0.5)
    cell_id = cell_subdf[1,:cell_id]
    for (i,p) in enumerate(levels(cell_subdf[!,:period]))
        idxs = cell_subdf.period .>= p
        # if all(cell_subdf[idxs,key] .> thresh)
        if (sum(cell_subdf[idxs,key] .> thresh)/length(cell_subdf[idxs,key])) > frac
            return i
        end
    end
    # hack so background is black with 16 stimulation conditions
    return -1
end
# combine(c->size(c,1),
#     groupby(avg_stim_per_period_df, :cell_id))

# thresh = 0.15
# thresh = 0.1 # kalman
# thresh = 0.0005
thresh = 0.001
# thresh = 0.0001
# thresh = -0.001
# df = avg_stim_per_period_df
df = avg_stim_df
stimthresh_df = combine(c->stim_threshold(c,:df_f, thresh),
    groupby(df, :cell_id))
# thresh = 1.00
# stimthresh_df = combine(c->stim_threshold(c,:df_f_mean, thresh),
#     groupby(avg_stim_df, :cell_id))

n1 = sum(stimthresh_df.x1 .== 1)
n2 = sum(stimthresh_df.x1 .== 2)
n3 = sum(stimthresh_df.x1 .== 3)
@show n1,n2,n3
@assert (n1>0) & (n2>0) & (n3>0)
combine(x->size(x,1), groupby(df, :period))
combine(groupby(df, :period), :df_f => mean)


idxs = stimthresh_df.x1 .>= 0
@show sum(idxs)
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

# need to sort by stim number 
lookup_dict = Dict(eachrow(recruited_stimthresh_df[:,[:cell_id, :x1]]))
recruited_labels = label_by_cellid(recruited_cell_masks, 
recruited_ids,
    tseriesH, tseriesW, tseriesZ,
    c->getindex(lookup_dict ,c));

@show sum(recruited_labels .> 0) 
## no need to run more than once!
error("don't run twice")
tseries_subset = tseries[:,:,:,1:20:tseriesT];
maxTproj = mapslices(x->quantile(x[:],0.90), tseries_subset,dims=4);
shiftedmaxTproj = shift_every_other_row(maxTproj,-3);
shiftedmaxTproj = shiftedmaxTproj[:,:,:,1];
##
shifted_rl = shift_every_other_row(recruited_labels,-3);
early_mid_late_im = RGB.(zeros(UInt8, size(shifted_rl)[1:2]...))
for i in 1:3
    is_recruited = shifted_rl.==i
    new_recruits = maximum(is_recruited, dims=3)[:,:,1]
    idxs = findall(new_recruits .> 0)
    if i == 1
        # Early is red
        r = RGB(0.8,0,0)
    elseif i == 2
        # mid is green
        r = RGB(0,0.8,0)
    else
        # late is yellow or blue
        # r = RGB(0.8,0.8,0)
        r = RGB(0,0,0.8)
    end
    early_mid_late_im[idxs] .= r
end

##
fig, ax = plt.subplots(1,dpi=600)
p = nothing
struct_im = maximum(shiftedmaxTproj, dims=3)[:,:,1]
# struct_im = struct_im[175:310,100:235]
struct_im = L.adjust_histogram(struct_im, GammaCorrection(0.5))
struct_im = RGB.(imadjustintensity(struct_im))
# eml_im = early_mid_late_im[175:310,100:235]
eml_im = early_mid_late_im
im = mapc.((x,y)->maximum([x,y]), eml_im, struct_im)
im = convert(Matrix{RGB{N0f16}}, im)
##
impath = joinpath(plot_dir,"early-mid-late_maxZ.png")
save(joinpath(plot_dir,"early-mid-late_maxZ.png"), im)
##
@pun trial_average = recording;
imap = influence_map(trial_average, window_len);
##
fig, ax = plt.subplots(1,dpi=600)
cmax = 0.1
cmin = -0.1
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
im = shift_every_other_row(imap[:,:,:,16],-3)
im = maximum(im, dims=3)[:,:,1]
cim = ax.imshow(im[175:310,100:235], cmap="RdBu_r",
            norm=cnorm)
cbar = plt.colorbar(cim)
ax.set_xticks([])
ax.set_xticks([], minor=true)
ax.set_yticks([])
ax.set_yticks([], minor=true)
fig

##

imshow(trial_average)