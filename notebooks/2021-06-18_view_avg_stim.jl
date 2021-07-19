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
plt = PyPlot
matplotlib = plt.matplotlib

##
ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
else
    # tseriesRootDir = "/data/dlab/b115"
    tseriesRootDir = "/scratch/b115"
    tseriesRootDir = "/mnt/deissero/users/tyler/b115"
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
##
avg_stim_h5_path = "$tseriesRootDir/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-titration-192trial-062_avgStim_lstm.h5"
avgStim = h5read(avg_stim_h5_path, "/block1");
##
imshow(avgStim)