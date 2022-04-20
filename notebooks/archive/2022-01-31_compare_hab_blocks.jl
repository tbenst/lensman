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