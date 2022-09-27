ENV["DISPLAY"] = "localhost:13.0"
##
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
using ImageView
##
# base_dir = "/data/dlab/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1"
base_dir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1"
tif_dirs = [
    "$(base_dir)/TSeries-resting-pre-038",
    "$(base_dir)/TSeries-resting-mid-040",
    "$(base_dir)/TSeries-resting-post-042",
    "$(base_dir)/TSeries-lrhab_raphe_40trial-039",
]

##
pre = avgTseries(tif_dirs[1]);
mid = avgTseries(tif_dirs[2]);
post = avgTseries(tif_dirs[3]);
stim1 = avgTseries(tif_dirs[4]);
##
imshow(stim1)