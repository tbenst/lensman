## compare 
# ENV["DISPLAY"] = "localhost:12.0"
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
# using ImageView
import Gadfly
using Unitful: μm, m, s, mW
import Base.Threads.@threads
import PyPlot
import Base.size
plt = PyPlot
matplotlib = plt.matplotlib

##
ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
else
    # tseriesRootDir = "/data/dlab/b115"
    tseriesRootDir = "/scratch/b115"
    # tseriesRootDir = "/mnt/deissero/users/tyler/b115"
end

# newer
if ON_SHERLOCK
    slmDir = "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files"
else
    slmDir = "/mnt/deissero/users/tyler/b115/SLM_files"
    # slmDir = "/mnt/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment"
    # slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
    # older
    # slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"
    # slmDir = "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/"
end

# tseriesDir = "$tseriesRootDir/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061" # ok
# tseriesDir = "$tseriesRootDir/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062" # ok
# tseriesDir = "$tseriesRootDir/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062" # ok
# tseriesDir = "$tseriesRootDir/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069" # ok
# tseriesDir = "$tseriesRootDir/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070" # ok
tseriesDir = "$tseriesRootDir/2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-cstoner-n64-b2-r8-077" # wrong number of z-planes...?

H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tseriesDir); # 15 z-planes...?
H, W, Z, T

##
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
Z = size(etlVals,1) # 11 z-planes...?

##
files = readdir(tseriesDir);
##
start = 999
files[start:(start+14)] # should be 10 planes!!!