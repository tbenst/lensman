ENV["DISPLAY"] = "localhost:10"
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
slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"

addprocs(36)
@everywhere begin
    import Pkg
    Pkg.activate(".")
    using Lensman
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
# to fix, run `mount -o remount,size=100G /dev/shm`
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> Lensman._readTseriesTiffDir(S, framePlane2tiffPath, Z, T));
##
voltageFile = glob("*VoltageRecording*.csv", dataDir)[1]
voltages = CSV.File(open(read, voltageFile)) |> DataFrame
rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)]);

##
# ! means get reference for df, not copy
plot(voltages[!,"Time(ms)"][1:10000], voltages[!,"frame starts"][1:10000])
frameStarted = diff(voltages[!,"frame starts"] .> std(voltages[!,"frame starts"]).*2)
frameStartIdx = findall(frameStarted .== 1) .+ 1;

##
@warn "hardcodes 10vol stim (2s @ 5Hz)"
stimStartIdx = findTTLStarts(voltages[!,"respir"])[1:10:end-9]
stimEndIdx = findTTLEnds(voltages[!,"respir"])[10:10:end];
##
stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimStartIdx]
stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimEndIdx]
# single stim example
imshow(imadjustintensity(tseries[:,:,:,stimStartFrameIdx[1]-1:stimEndFrameIdx[1]+1]))

## read slm stim files
dataFolders = splitpath(dataDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
open(xmlPath, "r") do io
    zseries_xml = read(io, String)
    global tseries_xml = xp_parse(zseries_xml)
end;
# warning: find has different syntax from xpath!!!
expDate = LibExpat.find(tseries_xml, "/PVScan[1]{date}")
expDate = DateTime(expDate, "m/d/y I:M:S p")
# much faster by prepending /PVScan/PVStateShard[1] to search first frame only
framePeriod = parse(Float64,tseries_xml[xpath"/PVScan/PVStateShard[1]/PVStateValue[@key='framePeriod']/@value"][1])
frameRate = 1/framePeriod
##
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))

"""Parse Sean's filename to time
trialOrder_2020_11_02___21h03m29sr.txt -> Time(21,3,29)"""
function getTimeFromFilename(fn)
    re = r".*_(?<year>\d+)_(?<month>\d+)_(?<day>\d+)___(?<hour>\d+)h(?<min>\d+)m(?<sec>\d+)sr?.txt$"
    gs = match(re, fn)
    dt = parse.(Int,(gs["year"], gs["month"], gs["day"],
                     gs["hour"], gs["min"], gs["sec"]))
    DateTime(dt...)
end

# only times that are within subfolders were used by a "Run"
# EXAMPLE (/mnt/deissero/users/tyler/b115/SLM_files/02-Nov-2020])
# e.g. after generating holograms, on first run (21h09m31s),
# a subfolder is created like 2020_11_2_21_6_50_, and Sean's code 
# copies trialOrder_2020_11_02___21h09m31s.txt to it.
#
# If we "Stop" and "Run" again without regenerating (but possibly
# changing trial order), no new folder is made BUT the new files for
# the RUN are copied to the same subfolder
# (trialOrder_2020_11_02___22h04m47s.txt)

slmRunDirs = joinpath.(slmExpDir,
    filter(x -> isdir(joinpath(slmExpDir, x)), readdir(slmExpDir)))

trialOrderTxts = vcat(glob.("trialOrder_*s.txt", slmRunDirs)...)
tempDts = getTimeFromFilename.(trialOrderTxts)

# sort by datetime
trialOrderTxts = trialOrderTxts[sortperm(tempDts)]
availableDateTimes = getTimeFromFilename.(trialOrderTxts)

# slmDate should always be younger than expDate
dtIdx = searchsortedfirst(availableDateTimes, expDate)-1
slmDateTime = availableDateTimes[dtIdx]
timeDiff = Dates.Second(expDate - slmDateTime)
println("started SLM $(timeDiff) seconds before imaging")
if (timeDiff > Second(120)) | (timeDiff < Second(0))
    @warn "timeDiff looks funky...verify we have correct trial"
end
trialOrderTxt = trialOrderTxts[dtIdx]

trialOrder = CSV.File(trialOrderTxt,delim="\t",header=true) |> Tables.matrix
@assert size(trialOrder,1)==1
trialOrder = trialOrder[1,:];
nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
@assert nTrials == size(stimStartIdx) # check TTL start-times match 

## avg stim effec
pre = Int(ceil(1*frameRate))
post = Int(ceil(1*frameRate))
trialTime = minimum(stimEndFrameIdx.-stimStartFrameIdx) + pre + post + 1
# HWZCT
avgStim = zeros(size(tseries)[1:3]..., nStimuli, trialTime)
##
for (i,(start, stop)) in enumerate(zip(stimStartFrameIdx, stimEndFrameIdx))
    trialType = trialOrder[i]::Int
    for (t,idx) in enumerate(start-pre:stop+post)
        avgStim[:,:,:,trialType,t] .+= tseries[:,:,:,idx]
    end
end
avgStim ./= nTrialsPerStimulus;
##
imshow(avgStim[:,:,:,1,:]) # trial type 1 (left hab)
imshow(avgStim[:,:,:,2,:]) # right hab
imshow(avgStim[:,:,:,3,:]) # raphe


## debug archive
# plot(voltages["Time(ms)"][1:110000], voltages["frame starts"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["secondary"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["winfluo"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["Blue"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["VR timestamps"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["green"][1:110000])
# plot(voltages["Time(ms)"][1:110000], voltages["LED"][1:110000])
# plot(voltages["Time(ms)"][1:120000], voltages["respir"][1:120000])