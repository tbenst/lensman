# ENV["DISPLAY"] = "localhost:12"
ENV["DISPLAY"] = ":0"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates, Distributed, SharedArrays
using Base.Iterators: peel

dataDir = "/scratch/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab_raphe_40trial-045/"

addprocs(18)
@everywhere using SharedArrays
@everywhere import ImageMagick


@everywhere function _readTseriesTiffDir(S::SharedArray, framePlane2tiffPath, Z, T)
    nProcs = size(procs(S),1)
    i = indexpids(S)
    myT = collect(i:nProcs:T)
    println("process $i will open $(size(myT,1)) tiffs")
    for t in myT
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            S[:,:,z,t] = ImageMagick.load(tp)
        end
    end
end
##
tifdir = dataDir
containsStr = "Ch3"
tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") &
    occursin(containsStr, x), readdir(tifdir)))
framePlane2tiffPath = Dict(getFramePlane(tp) => tp for tp in tiff_files)
framePlanes = hcat(collect.(keys(framePlane2tiffPath))...)
T = maximum(framePlanes[1,:])
Z = maximum(framePlanes[2,:])
tif0 = ImageMagick.load(tiff_files[1])
H, W = size(tif0)

# tseries = zeros(Normed{UInt16,16}, H, W, Z, T)
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> _readTseriesTiffDir(S, framePlane2tiffPath, Z, T))
# init = readTseriesTiffDir(framePlane2tiffPath, Z, T))

##
imshow(imadjustintensity(tseries))
# rawVolumes = r
##
tifPath = joinpath(dataDir, "TSeries-lrhab_raphe_40trial-045_Cycle03635_Ch3_000001.ome.tif")
tifPath = joinpath(dataDir, "TSeries-lrhab_raphe_40trial-045_Cycle03635_Ch3_000005.ome.tif")
tif = imadjustintensity(ImageMagick.load(tifPath))


##
tifdir = dataDir
contains = "Ch3"
tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") & occursin(contains, x),
    readdir(tifdir)))

##
function getFramePlane(tifPath)
    m = match(r".*Cycle(\d+)_Ch\d_(\d+).ome.tif", tifPath)
    parse(Int,m[1]), parse(Int,m[2])
end
frameNum, planeNum =getFramePlane(tifPath)
##
tif0 = ImageMagick.load(tiff_files[1])
H, W = size(tif0)

zseries = zeros(Normed{UInt16,16}, H, W, size(tiff_files, 1))
@threads for z in 1:size(tiff_files,1)
    zseries[:,:,z] = ImageMagick.load(tiff_files[z])
end
zseries


##
using Distributed
nprocs = 12
addprocs(nprocs)
@everywhere using SharedArrays
@everywhere import ImageMagick
##
H = 4; W = 4; Z = 2; T = 50
S = SharedArray{Normed{UInt16,16},4}(H, W, Z, T)
##
@distributed for t in 1:nprocs
    nProcs = size(procs(S),1)
    i = indexpids(S)
    myIdxs = collect(i:nProcs:T)
    println(myIdxs)
end
##
r