ENV["DISPLAY"] = "localhost:13"
##
# using AlgebraOfGraphics, CairoMakie
# aog = AlgebraOfGraphics
# Data = aog.data
# set_aog_theme!() #src

using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Dates, Arrow, StatsBase
using Gadfly
import Lensman: @pun, @assign
# using ImageView # sometimes hard to load over ssh
using AxisArrays
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
# init_workers(36)
Recordings
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



resources = Resources();
# wt_recordings = [Recordings[r](;resources...) for r in wt_tseries]
# rs_recordings = [Recordings[r](;resources...) for r in rs_tseries];

# annoyingly, can't run in parallel due to lazyTyH5 / HDF5 limits...
vol_rates = map(r-> r[:vol_rate], vcat(wt_recordings..., rs_recordings...))
proper_rate = isapprox.(vol_rates, 3; atol=0.1)
@assert all(proper_rate) vcat(wt_tseries..., rs_tseries...)[findall((~).(proper_rate))]
nTs = map(r-> r[:tseriesT], vcat(wt_recordings..., rs_recordings...))
@assert all(nTs .> 60) nTs
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

rs_ages = [est_age(r[:exp_date]) for r in rs_recordings]
wt_ages = [est_age(r[:exp_date]) for r in wt_recordings]
mean(rs_ages), mean(wt_ages)
##

wt_recordings[1][:uri]
wt_recordings[1][:tseries_dset]
# ramp_df = recording_to_df_df(wt_recordings[3], "wt")
ramp_df = vcat(collect(map(r->L.recording_to_df_df(r,"wt"), wt_recordings))...)
ramp_df = vcat(ramp_df, collect(map(r->L.recording_to_df_df(r,"rs"), rs_recordings))...)
open(joinpath(save_dir, "ramp_df.arrow"), "w") do io
    Arrow.write(io, ramp_df)
end;


##

the_count = combine(groupby(ramp_df, ["rel_time"]), :df_f => length)
@show the_count
##
plus_sem(x) = mean(x)+sem(x)
minus_sem(x) = mean(x)-sem(x)

df = ramp_df[ramp_df.fish .!= "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061", :]

meanDF = combine(groupby(df[df.rel_time .<= 58,:], ["rel_time", "genotype"]),
    "df_f" => mean, "df_f" => plus_sem, "df_f" => minus_sem, "rel_time" => r->r/3)
rename!(meanDF, "rel_time_function" => "time")
##
# Gadfly.plot(meanDF, x=:rel_time, ymin=:df_f_minus_sem, color=:genotype, ymax=:df_f_plus_sem, Geom.ribbon)
# Gadfly.plot(ramp_df, x=:rel_time, y=:df_f, color=:genotype, Geom.point)
# Gadfly.plot(df, x=:rel_time, y=:df_f, color=:genotype, Geom.smooth)
pmean = Gadfly.layer(meanDF, x=:time, y=:df_f_mean, color=:genotype, Geom.line)
pband = Gadfly.layer(meanDF, x=:time, ymin=:df_f_minus_sem, ymax=:df_f_plus_sem,color=:genotype, Geom.ribbon)
Gadfly.plot(pmean, pband)
# Gadfly.plot(ramp_df, x=:rel_time, y=:df_f, color=:fish, Geom.point, Guide.colorkey(labels=string.(collect(1:9))))
##
layers = smooth(span=0.2, degree=1)
p1 = data(meanDF) * (
    # mapping(:time, :df_f_mean, color=:genotype) * smooth(span=0.2, degree=1)
    # + mapping(:time, :df_f_minus_sem, :df_f_plus_sem, color=:genotype) * visual(Band)
    mapping(:time, :df_f_minus_sem, :df_f_plus_sem, color=:genotype) * visual(Band)
    # + mapping(:time, :df_f_minus_sem, color=:genotype)
    # + mapping(:time, :df_f_plus_sem, color=:genotype)
    # + mapping(:time, lower=:df_f_minus_sem, upper=:df_f_plus_sem, color=:genotype) * visual(LinesFill)
) |> draw
# data(ramp_df) * mapping(:time, :df_f, color=:genotype) * layers  |> draw
##
X = zeros(5,4,3,2)
Base.axes(X,2)

##
Xs = repeat(collect(1:100), 10)
Ys = sin.(Xs) + randn(1000)*0.2
df = DataFrame(x=Xs, y=Ys)
df = transform(groupby(df, ["x"]),
    "y" => mean, "y" => minimum, "y" => maximum)
data(df) * mapping(:x, :y_minimum, :y_maximum) * visual(Band) |> draw