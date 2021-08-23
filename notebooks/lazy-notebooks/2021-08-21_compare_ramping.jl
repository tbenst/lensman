ENV["DISPLAY"] = "localhost:13"
##
using AlgebraOfGraphics, CairoMakie
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Dates, Folds
import Lensman: @pun, @assign
using ImageView
using AxisArrays
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
aog = AlgebraOfGraphics
Data = aog.data
set_aog_theme!() #src
# init_workers(36)
Recordings
## find possible fish
data_dir = "/data/dlab/b115/"
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
    # "2021-03-16_rschrmine_h2b6s/fish3"
    # "2021-05-18_rsChRmine_h2b6s_6dpf/fish1"
    # "2021-05-18_rsChRmine_h2b6s_6dpf/fish3"
    "2021-05-18_rsChRmine_h2b6s_6dpf/fish5"
    "2021-06-01_rsChRmine_h2b6s/fish3"
    "2021-06-02_rsChRmine-h2b6s/fish2"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish2"
]

 wt = [
    # "2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2"
    # "2020-12-16_h2b6s-chrmine-kv2.1_6dpf/fish2"
    # "2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine"
    "2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish1_chrmine"
    "2021-02-16_h2b6s_wt-chrmine/fish3"
    # "2021-03-09_h2b6s/fish2"
    # "2021-03-16_h2b6s/fish1"
    # "2021-03-16_wt-chrmine_h2b6s/fish2"
    "2021-03-30_wt-chrmine_6dpf_h2b6s/fish1"
    "2021-04-13_wt-chrmine_6dpf_h2b6s/fish1"
    "2021-06-01_wt-chrmine_h2b6s/fish4"
]

function find_tseries(fish_dir; data_dir=data_dir)
    dirs = readdir(joinpath(data_dir,fish_dir))
    dirs = joinpath.(data_dir, fish_dir, dirs)
    tseries = filter(f->occursin("TSeries", f) && isdir(f), dirs)[end]
    tseries = splitpath(tseries)[end]
    if length(tseries) > 0
        return joinpath(fish_dir, tseries)
    end
end

wt_tseries = find_tseries.(wt)
rs_tseries = find_tseries.(rs)

resources = Resources();
wt_recordings = [Recordings[r](;resources...) for r in wt_tseries]
rs_recordings = [Recordings[r](;resources...) for r in rs_tseries];
wt_recordings[1][:uri]
wt_recordings[1][:tseries_dset]
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
##
# rs_ages = [est_age(r[:exp_date]) for r in rs_recordings]
# wt_ages = [est_age(r[:exp_date]) for r in wt_recordings]
# mean(rs_ages), mean(wt_ages)
##
@pun (tseries_dir, tyh5_path) = wt_recordings[2]
##
# wt_recordings[1] needs tyh5: /imaging/raw

# h5 = h5open(tyh5_path,"r")
lt = LazyTiff(tseries_dir)
if size(lt,3) != 10
    @warn "got $(size(lt,3)) z-planes"
end
##
Gray.(mean(imadjustintensity(lt[:,:,:,1:10]), dims=4))[:,:,1,1]

##
function tseries_to_df_df(tseries, genotype, fish; vol_rate = 3, f0_end=6,
     f1_end=72
)
    df = DataFrame()
    f0 = mean(tseries[:,:,:,1:f0_end])
    f = mean(tseries[:,:,:,1:f1_end], dims=[1,2,3])[1,1,1,:]
    rel_time = collect(1:length(f))
    time = rel_time./vol_rate
    df_f = @. (f-f0)/f0
    df = DataFrame(
        rel_time=rel_time,
        time= time,
        df_f= df_f,
        genotype=genotype,
        fish=fish
    )
end

function recording_to_df_df(recording, genotype)
    try
        @pun (tseries, vol_rate, uri) = recording
        f0_end = Int(round(vol_rate * 2))
        f1_end = Int(round(vol_rate * 20))
        fish = uri
        tseries_to_df_df(tseries, genotype, fish; vol_rate = vol_rate, f0_end=f0_end,
            f1_end=f1_end)
    catch
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println(stdout)
        end
        DataFrame()
    end
end

wt_recordings = [Recordings[r](;resources...) for r in wt_tseries]
rs_recordings = [Recordings[r](;resources...) for r in rs_tseries];

wt_recordings[1][:uri]
wt_recordings[1][:tseries_dset]
# ramp_df = recording_to_df_df(wt_recordings[1], "wt")
ramp_df = vcat(collect(Folds.map(r->recording_to_df_df(r,"wt"), wt_recordings)))
ramp_df = vcat(ramp_df, collect(Folds.map(r->recording_to_df_df(r,"rs"), rs_recordings)))

##
data(ramp_df) * mapping(:rel_time, :df_f)  |> draw
##
X = zeros(5,4,3,2)
Base.axes(X,2)