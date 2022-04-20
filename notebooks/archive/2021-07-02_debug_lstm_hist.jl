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
Geom = Gadfly.Geom

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


working_tseries_dir = "$tseriesRootDir/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"

# fixed with divide by 512
# not_working_tseries_dir = "$tseriesRootDir/2021-06-01_rsChRmine_h2b6s/fish3/TSeries-lrhab-118trial-060"

# too bright..?
not_working_tseries_dir = "$tseriesRootDir/2021-06-29_hsChRmine_6f_6dpf/fish2/TSeries-round3-lrhab-118trial-068"

##
function loadTseriesSamples(tifdir, containsStr::String="Ch3", num_frames=100)
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
    tseries = Array{UInt16}(undef, H, W, Z, num_frames)
    memory_size_bytes = prod(size(tseries)) * 2
    memory_size_gb = round(memory_size_bytes / 1024^3, digits=1)
    println("estimated memory usage: $memory_size_gb")
    p = Progress(T, 1, "Load Tseries: ")
    idxs = randperm(T)
    @threads for t in 1:num_frames
        tt = idxs[t]
        for z in 1:Z
            tp = framePlane2tiffPath[tt,z]
            tseries[:,:,z,t] .= reinterpret(UInt16, ImageMagick.load(tp))
        end
        next!(p)
    end
    tseries
end

workingtseries = loadTseriesSamples(working_tseries_dir,"Ch3",100)
badtseries = loadTseriesSamples(not_working_tseries_dir,"Ch3",100);
##
k = 10
nW = length(workingtseries[1:k:end])
nB = length(badtseries[1:k:end])

# Gadfly.plot(x=rand(100),
#             color=vcat(repeat(["blue"],50),repeat(["white"],50)),
#             Geom.histogram(bincount=50))

pw = Gadfly.plot(x=log.(workingtseries[1:k:end]),
            Geom.histogram(bincount=50),
            color=repeat(["working"],nW))
pb = Gadfly.plot(x=log.(badtseries[1:k:end]),
            Geom.histogram(bincount=50),
            color=repeat(["bad"],nW))

vstack(pw,pb)

