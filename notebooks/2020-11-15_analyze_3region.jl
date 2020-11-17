ENV["DISPLAY"] = "localhost:11"
# ENV["DISPLAY"] = ":0"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics
using Base.Iterators: peel

dataDir = "/scratch/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab_raphe_40trial-045/"

# where stim targets are defined
targetDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2"

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
# warning: this may exhaust shared memory, causing a bus error
# to fix, run `mount -o remount,size=90G /dev/shm`
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> _readTseriesTiffDir(S, framePlane2tiffPath, Z, T))
##
voltageFile = glob("*VoltageRecording*.csv", dataDir)[1]
voltages = CSV.File(open(read, voltageFile)) |> DataFrame
rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

##
plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
frameStarted = diff(voltages["frame starts"] .> std(voltages["frame starts"]).*2)
frameStartIdx = findall(frameStarted .== 1) .+ 1

##
@warn "hardcodes 10vol stim (2s @ 5Hz)"
stimStartIdx = findTTLStarts(voltages["respir"])[1:10:end-9]
stimEndIdx = findTTLEnds(voltages["respir"])[10:10:end]
##
stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimStartIdx]
stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimEndIdx]
imshow(imadjustintensity(tseries[:,:,:,stimStartFrameIdx[1]-1:stimEndFrameIdx[1]+1]))
## avg stim effect
imshow()



# plot(voltages["Time(ms)"][1:110000], voltages["frame starts"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["secondary"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["winfluo"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["Blue"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["VR timestamps"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["green"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["LED"][1:110000])
# plot(voltages["Time(ms)"][1:120000], voltages["respir"][1:120000])