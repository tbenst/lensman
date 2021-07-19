using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
using ImageView
import Gadfly
using Unitful: μm, m, s, mW
import Base.Threads.@threads
import PyPlot
import Base.size
plt = PyPlot
matplotlib = plt.matplotlib

ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
    slmDir = "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files"
else
    # tseriesRootDir = "/data/dlab/b115"
    tseriesRootDir = "/scratch/b115"
    # tseriesRootDir = "/mnt/deissero/users/tyler/b115"
    slmDir = "/mnt/deissero/users/tyler/b115/SLM_files"
end


tseriesDir = "$tseriesRootDir/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"


fish_dir = joinpath(splitpath(tseriesDir)[1:end-1]...)
plot_dir = joinpath(fish_dir, "plots-denoised")
if ~isdir(plot_dir)
    mkdir(plot_dir)
end

##
# tyh5Path = tseriesDir * ".ty.h5"
# tyh5Path = tseriesDir * "2021-06-21_6pm.ty.h5"
tyh5Path = tseriesDir * "skipme"
if isfile(tyh5Path)
    @info "using ty.h5 file"
    # tseries = read_tyh5(tyh5Path)
    # h5 = h5open(tyh5Path,"r")
    dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm"
    # dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden_init_from_pretrained-2021-06-21_6pm"
    # h5, tseries = lazy_read_tyh5(tyh5Path, dset);
    h5, tseries = read_tyh5(tyh5Path, dset);
    plotDir = joinpath(fishDir, "plots-denoised")
    avgStimStr = "_avgStim_lstm.h5"
else
    # error("stop")
    tseries = loadTseries(tseriesDir);
    plotDir = joinpath(fishDir, "plots")
    avgStimStr = "_avgStim.h5"
end;
##
function make_ellipse_mask(x,y,a,b,thesize=size(tseries)[1:2])
    mask = Gray.(zeros(Bool,thesize))
    draw!(mask, Ellipse(
        Point(x,y), a, b))
    reinterpret(Bool,mask)
end
rhab_mask = make_ellipse_mask(388,268,25,25)
lhab_mask = make_ellipse_mask(388,188,25,25)
hab_mask = rhab_mask .| lhab_mask;
hab_mask = repeat(hab_mask,1,1,10);
hab_mask = findall(channelview(hab_mask))
habTrace = extractTrace(tseries, hab_mask);
##
plot(habTrace)