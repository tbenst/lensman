##
ENV["DISPLAY"] = "localhost:13.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random
import Gadfly
using Unitful: μm, m, s

##
offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
zOffset = offset * 1e6

tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-047.ty.h5"
fishDir = joinpath(splitpath(tyh5Path)[1:end-1]...)
expName = replace(splitpath(tyh5Path)[end], ".ty.h5" => "")
tables = pyimport_conda("tables", "pytables")
tyh5 = tables.open_file(tyh5Path, "r")

# performance is thrashing, sometimes down to <100KB/s
# seems to be a filesystem issue

# pytables seems to read at ~35MB/s
tseries = tyh5.root.imaging.raw.read();
tseries = permutedims(tseries,(4,5,3,2,1));
tseries = tseries[:,:,:,1,:];
# HDF5.jl seems to read at ~20MB/s
# tseries = h5read(tyh5Path, "/imaging/raw");