ENV["DISPLAY"] = "localhost:13"
##
# using AlgebraOfGraphics, CairoMakie
# aog = AlgebraOfGraphics
# Data = aog.data
# set_aog_theme!() #src

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
Recordings;
## find possible fish
data_dir = "/data/dlab/b115/"
save_dir = "/data/dlab/b115/multi-fish/cross-stim"
# fd_path = "/nix/store/bpi6qsm9i7p9af8lk1snjjqi9wng0bpm-user-environment/bin/fd"
# cmd = `$fd_path -t d -d 3 TSeries $data_dir`

# nparent_dir = length(splitpath(data_dir))
# folders = splitpath.(split(read(cmd, String),"\n"))[1:end-1]
# ndirs = lastindex(folders[1])
# folders = [f[5:end] for f in folders]
# fish = map(p->joinpath(p...), (unique([f[1:2] for f in folders])))
# rs = filter(p->occursin("rs",p), fish)
# wt = filter(p->~(occursin("rs",p) | occursin("hs",p) | occursin("h33r",p)), fish)
# wt = filter(p->occursin("h2b6s",p), wt)
# rs[3:end]
# wt[end-7:end]
##
rs = [
    "2021-03-16_rschrmine_h2b6s/fish3"
    "2021-05-18_rsChRmine_h2b6s_6dpf/fish1"
    "2021-05-18_rsChRmine_h2b6s_6dpf/fish3"
    "2021-05-18_rsChRmine_h2b6s_6dpf/fish5"
    "2021-06-01_rsChRmine_h2b6s/fish3"
    "2021-06-02_rsChRmine-h2b6s/fish2"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish2"
]

 wt = [
    # "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1" # 2.5Hz
    # "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2" # 1.9Hz
    "2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1"
    "2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2"
    "2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2"
    "2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine"
    "2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish1_chrmine"
    # "2021-02-16_h2b6s_wt-chrmine/fish3" # 30Hz
    # "2021-03-09_h2b6s/fish2" #gcamp?
    # "2021-03-16_h2b6s/fish1" #gcamp?
    # "2021-03-16_wt-chrmine_h2b6s/fish2" # could specify TSeries_cross-stim_p125_10plane-211
    # "2021-03-30_wt-chrmine_6dpf_h2b6s/fish1" # 30Hz
    # "2021-04-13_wt-chrmine_6dpf_h2b6s/fish1" # 30Hz
    "2021-06-01_wt-chrmine_h2b6s/fish4"
]

gcamp_tseries = [
    # "2021-03-09_h2b6s/fish2/TSeries-32cell-8concurrent-10MHz-8rep-066" # 30Hz
    "2021-03-16_h2b6s/fish1/TSeries_cross-stim_p125_10plane-203" #gcamp?
    "2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish2_nochrmine/TSeries-1024cell-32concurrent-4power-048"
]

function find_tseries(fish_dir; data_dir=data_dir)
    dirs = readdir(joinpath(data_dir,fish_dir))
    dirs = joinpath.(data_dir, fish_dir, dirs)
    # tseries = filter(f->occursin("TSeries", f) && isdir(f), dirs)[end]
    tseries = filter(f->occursin("TSeries", f) && isdir(f), dirs)[1]
    tseries = splitpath(tseries)[end]
    if length(tseries) > 0
        return joinpath(fish_dir, tseries)
    end
end

wt_tseries = find_tseries.(wt)
rs_tseries = find_tseries.(rs)
function swap_tseries(tseries_strs, pair)
    idx = findall(tseries_strs .== pair[1])[1]
    tseries_strs[idx] = pair[2]
    tseries_strs
end

wt_tseries = swap_tseries(wt_tseries, "2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2/TSeries-1024cell-32concurrent-039"
    => "2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2/TSeries-1024cell-32concurrent-042")
rs_tseries = swap_tseries(rs_tseries, "2021-05-18_rsChRmine_h2b6s_6dpf/fish1/TSeries-lrhab-control-91trial-4Mhz-044"
    => "2021-05-18_rsChRmine_h2b6s_6dpf/fish1/TSeries-lrhab-control-91trial-4Mhz-047")
rs_tseries = swap_tseries(rs_tseries, "2021-03-16_rschrmine_h2b6s/fish3/TSeries_64cell_8concurrent_2power_8rep-407"
    => "2021-03-16_rschrmine_h2b6s/fish3/TSeries_cross-stim_p125_10plane-313")
rs_tseries = swap_tseries(rs_tseries, "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-cstoner-n64-b2-r8-077"
    => "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069")

# exclusion criteria:
# - want at least 2 min post previous recording
# - age of fish should be 5-7dpf
# - baseline expression should be matched between cohorts

wt_tseries = vcat(wt_tseries...,
    "2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2/TSeries-1024cell-32concurrent-042",
    "2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-resting-pre-037",
    "2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_resting_pre-042",
    # "2020-10-27_elavl3-chrmine-Kv2.1_h2b6s_7dpf/fish1/TSeries-resting-pre-037", # 1.6Hz
    "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-resting-pre-029",
    "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1/TSeries-resting-pre-019", # quickly followed another?
    "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-resting-pre-029",
    # "2020-10-25_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-post-017"
)

resources = Resources();
wt_recordings = [Recordings[r](;resources...) for r in wt_tseries]
rs_recordings = [Recordings[r](;resources...) for r in rs_tseries];
gcamp_recordings = [Recordings[r](;resources...) for r in gcamp_tseries];


recordings = vcat(wt_recordings..., rs_recordings..., gcamp_recordings...)
recording_genotypes = vcat(
    repeat(["wt"], length(wt_recordings)),
    repeat(["rs"], length(rs_recordings)),
    repeat(["gcamp"], length(gcamp_recordings))
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
##
t_vol_rates = Folds.map(r-> r[:vol_rate], threads_ok_recordings)
nt_vol_rates = map(r-> r[:vol_rate], nothreads_recordings)
vol_rates = vcat(t_vol_rates..., nt_vol_rates...)
proper_rate = isapprox.(vol_rates, 3; atol=0.1)
@assert all(proper_rate) uris[findall((~).(proper_rate))]
##
t_TseriesTs = Folds.map(r-> r[:tseriesT], threads_ok_recordings)
nt_TseriesTs = map(r-> r[:tseriesT], nothreads_recordings)
TseriesTs = vcat(t_TseriesTs..., nt_TseriesTs...)
@assert all(TseriesTs .> 60) TseriesTs
##
##
# recording = Recordings[joinpath(data_dir,
#     "2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-resting-pre-029")](;
#     resources...)
# @pun (vol_rate, tseriesT) = recording
##
function est_age(date)
    day = dayname(date)
    val = 0
    if day == "Monday"
        val += 5
    elseif day == "Tuesday"
        val += 6
    elseif day == "Wednesday"
        val += 7
    elseif day == "Thursday"
        val += 8
    elseif day == "Friday"
        val += 9
    elseif day == "Saturday"
        val += 10
    elseif day == "Sunday"
        val += 4
    end
    h = hour(date) - 9 # fish typically spawn at 9am
    round(val + h/24,digits=2)
end

ages = [est_age(r[:exp_date]) for r in recordings]
# mean(rs_ages), mean(wt_ages)
##

ramp_df1 = vcat(collect(Folds.map(x->L.recording_to_df_df(x[1],x[2]),
    zip(threads_ok_recordings, thread_ok_genotypes)))...)
ramp_df2 = vcat(collect(map(x->L.recording_to_df_df(x[1], x[2]),
    zip(nothreads_recordings, nothread_genotypes)))...)
ramp_df = vcat(ramp_df1, ramp_df2)
# ramp_df = recording_to_df_df(wt_recordings[3], "wt")
# ramp_df = vcat(collect(map(r->L.recording_to_df_df(r,"wt"), wt_recordings))...)
# ramp_df = vcat(ramp_df, collect(map(r->L.recording_to_df_df(r,"rs"), rs_recordings))...)
open(joinpath(save_dir, "ramp_df_f0=1_eps=1.arrow"), "w") do io
    Arrow.write(io, ramp_df)
end;


##

the_count = combine(groupby(ramp_df, ["rel_time"]), :df_f => length)
@show the_count
##
plus_sem(x) = mean(x)+sem(x)
minus_sem(x) = mean(x)-sem(x)
# plus_sem(x) = quantile(x,0.75)
# minus_sem(x) = quantile(x,0.25)
df = ramp_df
##
@show df[df.rel_time .== 1, [:f, :genotype, :fish]]

fs = df[df.rel_time .== 1, :f]
fishes = df[df.rel_time .== 1, :fish]
outlier_upper_thresh = mean(fs) + 2*std(fs)
outlier_lower_thresh = mean(fs) - 2*std(fs)
outlier_fish = fishes[(fs .> outlier_upper_thresh) .| (fs .< outlier_lower_thresh)]
##
# we keep first 58 to avoid stim artifact / account for vol_rate differences
df = ramp_df[map(f-> ~in(f, outlier_fish), ramp_df.fish), :]
df = df[df.genotype .!= "gcamp", :] # N is too low to 
meanDF =combine(groupby(df[df.rel_time .<= 58,:], ["rel_time", "genotype"]),
    "df_f" => mean, "df_f" => plus_sem, "df_f" => minus_sem)
meanDF[:,:time] = meanDF[:,:rel_time]/mean(vol_rates)
combine(groupby(df, :genotype), :fish => f->length(unique(f)))

## Plot kinetics of cross-stim
# Gadfly.plot(meanDF, x=:rel_time, ymin=:df_f_minus_sem, color=:genotype, ymax=:df_f_plus_sem, Geom.ribbon)
# Gadfly.plot(df, x=:rel_time, y=:df_f, color=:genotype, Geom.smooth)

W = 183mm / 2
H = W * 0.8
set_default_graphic_size(W, H)
# ppoints = Gadfly.layer(ramp_df, x=:time, y=:df_f, xgroup=:genotype, color=:genotype, Geom.point)
pmean = Gadfly.layer(meanDF, x=:time, y=:df_f_mean, xgroup=:genotype, color=:genotype,
    Geom.line)
pband = Gadfly.layer(meanDF, x=:time, ymin=:df_f_minus_sem, ymax=:df_f_plus_sem,
    xgroup=:genotype, color=:genotype, Geom.ribbon)
# Gadfly.plot(Geom.subplot_grid(pmean, pband, ppoints))

p_cross = Gadfly.plot(Geom.subplot_grid(pmean, pband), Guide.xlabel("time (s)"),
    Guide.ylabel("Δf/f"))
# Gadfly.plot(ramp_df, x=:rel_time, y=:df_f, color=:fish, Geom.point, Guide.colorkey(labels=string.(collect(1:9))))
# Gadfly.plot(df[df.rel_time .== 58, :], x=:genotype, y=:df_f, color=:genotype, Geom.boxplot)
svg_fn = joinpath(save_dir, "16fish-cross-stim_rs-vs-wt.svg")
fig = draw(SVG(svg_fn), p_cross)
L.inkscape_make_png(svg_fn, H, W)
fig

##
df_last_seconds = df[(df.rel_time .> 51 .& df.rel_time .<= 57),:]
# df_last_seconds = df[(df.rel_time .> 56 .& df.rel_time .<= 57),:]
# df_last_seconds = df[df.rel_time .== 54 ,:]
x = df_last_seconds[(df_last_seconds.genotype .== "wt"),:df_f]
y = df_last_seconds[(df_last_seconds.genotype .== "rs"),:df_f]
HypothesisTests.MannWhitneyUTest(x,y)
##
meanDF_last_secs =combine(groupby(df_last_seconds, "genotype"),
    "df_f" => mean, "df_f" => plus_sem, "df_f" => minus_sem)
meanDF_last_secs[meanDF_last_secs.genotype .== "wt",:genotype] .= "ChRmine"
meanDF_last_secs[meanDF_last_secs.genotype .== "rs",:genotype] .= "rsChRmine"

Gadfly.plot(meanDF_last_secs, x=:genotype, y=:df_f_mean, ymin=:df_f_minus_sem, ymax=:df_f_plus_sem,
    color=:genotype, Geom.point, Geom.errorbar,
    Guide.ylabel("Δf/f"))
##