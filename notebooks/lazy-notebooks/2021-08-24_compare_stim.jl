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
    "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-lrhab_raphe_stim-40trial-proper-zoffset-034",
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
open(joinpath(save_dir, "cells_df_f.arrow"), "w") do io
    Arrow.write(io, cells_df_f_combined)
end;
# TODO: need to add "fish" & "genotype" to each row. can count rows in the
# non-vcat version. Also, some stims are outside of the brain and those cells
# need to be excluded. can look at x,y location as a heuristic for each group
# of stims..? or just always select stim_group == 1 or 2?

##
recording = recordings[1];
getExpData(recording[:tseries_xml]) # etl_vals is nothing???
@pun (
    cells, tseriesH, tseriesW, tseriesZ, target_size_px, vol_rate,
    tseries
) = recording;

cell_masks = constructROImasks(cells, tseriesH, tseriesW, tseriesZ, target_size_px)
##
cells_df_f_win_secs = 2
cells_df_f_delay = 0
cells_df_f_winsize = Int(ceil(cells_df_f_win_secs*vol_rate));
# TODO check if this worked..?
# was 16minutes without threads...
# with threads it's ~20min even though we are slamming CPU...?
cells_df_f = L.add_df_f_to_cells(tseries, cells, cell_masks,
    winSize=cells_df_f_winsize, delay=cells_df_f_delay, threads=false);
##
cells_df_f # TODO: why is stimFeq an array..?
# also TODO why are some df_f values crazy outlier
cells_df_f[:,7:10]
##
cells_df_f[:,:df_f_symlog] .= L.symlog.(cells_df_f[:,:df_f])
##
plot(cells_df_f, y=:df_f_symlog, Geom.violin)
##
mean_cells = combine(groupby(cells_df_f,:cellID), "df_f" => mean)
mean_cells[:,:df_f_symlog] .= L.symlog.(mean_cells[:,:df_f_mean])

plot(mean_cells, y=:df_f_symlog, Geom.violin)