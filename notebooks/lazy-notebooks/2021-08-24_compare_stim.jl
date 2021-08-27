ENV["DISPLAY"] = "localhost:13"
##
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Dates, Arrow, StatsBase, Folds, HypothesisTests
using Gadfly
import Lensman: @pun, @assign
# using ImageView # sometimes hard to load over ssh
using AxisArrays
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
# init_workers(36)
## choose fish
data_dir = "/data/dlab/b115/"
save_dir = "/data/dlab/b115/multi-fish/cross-stim"


wt_tseries = vcat(
    # 1Mhz
    # VoltageRecording file returns no stim times..?
    # "2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-039",
    "2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab-raphe-20trial-part2-047",
    "2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-lrhab_raphe_stim-40trial-038",
    "2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_lrhab_raphe_40trial-044",
    # making cells_df_f array fails..?
    # "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-lrhab_raphe_stim-40trial-proper-zoffset-034",
    "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1/TSeries-lrhab_raphe_40trial-023",
    # was this 4Mhz? also, started slm 1100 seconds before...suspicious...
    # see /mnt/deissero/users/tyler/b115/SLM_files/18-May-2021
    "2021-04-13_wt-chrmine_6dpf_h2b6s/fish1/TSeries-4region-lrhab-raphe-control-129trial-145",

    "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
)

rs_tseries = vcat(
    "2021-05-18_rsChRmine_h2b6s_6dpf/fish1/TSeries-lrhab-control-91trial-4Mhz-047",
    # forgot VoltageRecording...?
    # "2021-05-18_rsChRmine_h2b6s_6dpf/fish3/TSeries-lrhab-control-91trial-4Mhz-057",
    "2021-05-18_rsChRmine_h2b6s_6dpf/fish5/TSeries-lrhab-control-91trial-4Mhz-045",
    "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-lrhab-118trial-060",
    "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061",
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061",
    # for some reason this also fails for df_f...?
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
)
gcamp_tseries = vcat(
    "2021-01-25_gcamp6f_6dpf/fish1/TSeries-1024cell-32concurrent-4power-044",
    "2021-02-09_gcamp6f_7dpf/fish2/TSeries-lrhab-raphe-balanced-transitions-137",
    "2021-02-23_h2b6s/fish1_pt2/TSeries-128cell-4concurrent-3power-049",
    "2021-03-09_h2b6s/fish2/TSeries-32cell-8concurrent-10MHz-8rep-066",
    "2021-03-16_h2b6s/fish1/TSeries_64cell_8concurrent_2power_8rep-207"
)

resources = Resources();
wt_recordings = [Recordings[r](;resources...) for r in wt_tseries]
rs_recordings = [Recordings[r](;resources...) for r in rs_tseries];
# gcamp_recordings = [Recordings[r](;resources...) for r in gcamp_tseries];


recordings = vcat(wt_recordings..., rs_recordings...)
recording_genotypes = vcat(
    repeat(["wt"], length(wt_recordings)),
    repeat(["rs"], length(rs_recordings)),
)
threads_ok = [isnothing(r[:tseries_dset]) for r in recordings]
threads_ok_recordings = recordings[threads_ok]
# annoyingly, can't run in parallel due to HDF5 limits...
nothreads_recordings = recordings[(~).(threads_ok)];
recordings = vcat(threads_ok_recordings..., nothreads_recordings...);
thread_ok_genotypes = recording_genotypes[threads_ok]
nothread_genotypes = recording_genotypes[(~).(threads_ok)]
recording_genotypes = vcat(thread_ok_genotypes..., nothread_genotypes...)
uris = [r[:uri] for r in recordings]

# df = recordings[1][:cells_df_f];
"done"
## make sure similar recording settings
vol_rates = L.thread_safe_map(r-> r[:vol_rate],
    threads_ok_recordings, nothreads_recordings)
proper_rate = isapprox.(vol_rates, 3; atol=0.1)
@assert all(proper_rate) uris[findall((~).(proper_rate))]
##
TseriesTs = L.thread_safe_map(r-> r[:tseriesT],
    threads_ok_recordings, nothreads_recordings)
@assert all(TseriesTs .> 60) TseriesTs
## quick check one recording
# recording = Recordings[joinpath(data_dir,
#     "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-resting-pre-029")](;
#     resources...)
# @pun (vol_rate, tseriesT) = recording
##
ages = [L.est_age(r[:exp_date]) for r in recordings]
@assert all(ages .>= 5)
@assert all(ages .< 8)
##
function get_stim_start(r)
    try
        r[:stim_start_idx]
    catch
        []
    end
end

stim_start_idxs = L.thread_safe_map(get_stim_start,
# stim_start_idxs = L.thread_safe_map(r->r[:stim_start_idx],
    threads_ok_recordings, nothreads_recordings);
uris[length.(stim_start_idxs) .== 0]
##
cells_df_f = L.thread_safe_map(r-> r[:cells_df_f],
    threads_ok_recordings, nothreads_recordings; noerror=true);
# for some reason, "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-lrhab_raphe_stim-40trial-proper-zoffset-034"
# fails
##
for c in 1:length(cells_df_f)
    if isnothing(cells_df_f[c])
        continue
    end
    cells_df_f[c][!, :genotype] .= recording_genotypes[c]
    cells_df_f[c][!, :uri] .= uris[c]
end
##
good_idxs = (~).(isnothing.(cells_df_f))
cells_df_f_combined = vcat(cells_df_f[good_idxs]...)
open(joinpath(save_dir, "cells_df_f_f0_nokalman.arrow"), "w") do io
    Arrow.write(io, cells_df_f_combined)
end;
# Also, some stims are outside of the brain and those cells
# need to be excluded. can look at x,y location as a heuristic for each group
# of stims..? or just always select stim_group == 1 or 2?
##
df = cells_df_f_combined[cells_df_f_combined.stimNum .<= 2, :]
##
tenpercentile(x) = quantile(x, 0.1)
median_df = combine(groupby(df, [:uri, :cellID, :genotype]),
    "df_f" => median, "f" => median, "f0" => tenpercentile)
@. median_df[:,:perc_df_f] = (median_df[:,:f_median] - median_df[:,:f0_tenpercentile]) / (median_df[:,:f0_tenpercentile]+10)

# this does the best job of denoising
sum_df = combine(groupby(df, [:uri, :cellID, :genotype]),
    "f" => sum, "f0" => sum)
region_df = combine(groupby(df, [:uri, :stimGroup, :genotype]),
    "f" => sum, "f0" => sum)
@. region_df[:,:region_df_f] = (region_df[:,:f_sum] - region_df[:,:f0_sum]) / (region_df[:,:f0_sum]+10)
region_df_err = combine(groupby(region_df, [:genotype]),
    "region_df_f" => mean,
    "region_df_f" => sem, "region_df_f" => L.plus_sem, "region_df_f" => L.minus_sem,
    "region_df_f" => std, "region_df_f" => L.plus_std, "region_df_f" => L.minus_std)
@show region_df_err[:,[:region_df_f_mean, :region_df_f_sem]]
n_wt_trials = size(region_df[region_df[:,:genotype].=="wt",:], 1)
n_rs_trials = size(region_df[region_df[:,:genotype].=="rs",:], 1)
n_rs_fish = size(unique(region_df[region_df[:,:genotype].=="rs",:uri]), 1)
n_wt_fish = size(unique(region_df[region_df[:,:genotype].=="wt",:uri]), 1)
@show n_wt_trials, n_rs_trials, n_rs_fish, n_wt_fish

##
# dff = median_df[median_df.df_f_median .>= 0.2, :]
dff = median_df
# plot(dff, x=:genotype, y=:df_f_median, Geom.boxplot,
#     Coord.cartesian(;ymin=-1, ymax=7))
# plot(sum_df, color=:genotype, x=:sum_df_f, Geom.density,)
# plot(sum_df, color=:genotype, x=:sum_df_f, Geom.density,)
# plot(region_df, color=:genotype, x=:region_df_f, Geom.density,)
plot(region_df, color=:genotype, x=:genotype, y=:region_df_f, Geom.boxplot,
    Coord.cartesian(;ymin=-1, ymax=5))
plot(region_df_err, x=:genotype, y=:region_df_f_mean,
    # ymin=:region_df_f_minus_sem, ymax=:region_df_f_plus_sem,
    ymin=:region_df_f_minus_std, ymax=:region_df_f_plus_std,
    color=:genotype, Geom.point, Geom.errorbar,
    Coord.cartesian(;ymin=-1, ymax=2))
##,
med_per_fish = combine(groupby(median_df, [:uri, :genotype]),
    "df_f_median" => median, "f_median" => median, "f0_median" => median)[:,2:end]
##
f = med_per_fish.f_median_median
f0 = med_per_fish.f0_median_median

@. (f-f0)/f
# .83 wt; .36 rs
##
# 4 Mhz
idxs = (
        median_df.uri .== "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
    ) .| (
        median_df.uri .== "2021-04-13_wt-chrmine_6dpf_h2b6s/fish1/TSeries-4region-lrhab-raphe-control-129trial-145"
    )

mean(median_df[idxs,:df_f_median])
# .27
##
mean_cells = combine(groupby(cells_df_f,:cellID), "df_f" => mean)
mean_cells[:,:df_f_symlog] .= L.symlog.(mean_cells[:,:df_f_mean])

plot(mean_cells, y=:df_f_symlog, Geom.violin)