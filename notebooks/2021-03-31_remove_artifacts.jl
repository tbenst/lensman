using FileIO, NRRD, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, Images, #ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, HDF5, Arrow, ProgressMeter,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman
import Base.Threads.@threads
using Gadfly
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
##
tif_dir = "/data/dlab/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-lrhab_raphe_stim-40trial-038"
##
function loadTseriesBOT(tifdir, containsStr="Ch3")
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
    tseries_bot = Array{Float64}(undef, Z, T)
    p = Progress(T, 1, "Load Tseries: ")
    T = 1000
    @threads for t in 1:T
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            tseries_bot[z,t] = mean(reinterpret(UInt16, ImageMagick.load(tp)))
        end
        next!(p)
    end
    tseries_bot
end

function loadTseriesSmall(tifdir, containsStr="Ch3", maxT=1000)
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
    T = 1000
    tseries = Array{UInt16}(undef, H, W, Z, T)
    memory_size_bytes = prod(size(tseries)) * 2
    memory_size_gb = round(memory_size_bytes / 1024^3, digits=1)
    println("estimated memory usage: $memory_size_gb")
    p = Progress(T, 1, "Load Tseries: ")

    @threads for t in 1:T
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            tseries[:,:,z,t] .= reinterpret(UInt16, ImageMagick.load(tp))
        end
        next!(p)
    end
    tseries
end
##
# loadTseriesSmall
green_bot = loadTseriesBOT(tif_dir);
red_bot = loadTseriesBOT(tif_dir, "Ch2");