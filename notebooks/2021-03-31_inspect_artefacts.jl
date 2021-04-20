# for looking at output of chris's scripts...
# WIP since output not working yet on my data :/

using FileIO, NRRD, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, Images, #ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, HDF5, Arrow,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman
using Gadfly
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
##

data_dir = "/data/dlab/b115/process-output/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-resting-pre-038"

data_dir*"hdf5"
data_dir*"output"
##
h5open(joinpath(data_dir,"output", "artefact.h5"), "r") do h5
    keys(h5)
    # @show size(h5["data"])
    # @show size(h5["stim_start"])
    # @show size(h5["stim_stop"])
    h5["stim_stop"]
end
##
tseries = loadTseries("/data/dlab/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-039");
