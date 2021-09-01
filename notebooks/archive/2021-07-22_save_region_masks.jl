using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions, NRRD, NIfTI, ProgressMeter
using ImageView
import Gadfly
using Unitful: μm, m, s, mW
import Base.Threads.@threads
import PyPlot
import Base.size
plt = PyPlot
matplotlib = plt.matplotlib
sio = pyimport("scipy.io")

ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
else
    # tseriesRootDir = "/scratch/b115"
    tseriesRootDir = "/mnt/deissero/users/tyler/b115/"
    zbrain_dir = "/mnt/deissero/users/tyler/zbrain"
end

fishDir = joinpath(tseriesRootDir,
    "2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1")

##
fixedname = joinpath(fishDir, "zseries.nrrd")
@assert isfile(fixedname)
zseries = load(fixedname);
##

affine_transform_path = glob("*GenericAffine.mat", fishDir)[1]

affine_transform = sio.loadmat(affine_transform_path)
masks = matread("$zbrain_dir/MaskDatabase.mat")
raphe_mask = getMaskNameIncludes(masks, "Raphe");
##
# make raphe binary again
raphe_mask = antsApplyTransforms(fixedname, raphe_mask, affine_transform_path) .> 0

##
function registeredmask(masks, idx, fixedname, affine_transform_path)
    name = mask_names[idx]
    mask = read_mask(masks, idx)
    mask = AxisArray(Float32.(mask), (:y, :x, :z), zbrain_units)
    antsApplyTransforms(fixedname, mask, affine_transform_path) .> 0
end
##
mask_names = masks["MaskDatabaseNames"][1,:]
nMasks = length(mask_names)
# h5path = joinpath(fishDir, "masks.h5")
h5path = joinpath("/tmp", "masks.h5")
# rm(h5path)
@assert ~isfile(h5path)
h5 = Lensman.HDF5ThreadSafe.hdf5server(h5path, "w")

p = Progress(nMasks, 1, "Register:")
@threads for i in 1:nMasks
    name = replace(mask_names[i], "/" => "_") # or else a h5 path..
    registered = registeredmask(masks, i, fixedname, affine_transform_path)
    h5[name] = collect(registered) # convert BitArray to Array{Bool}
    next!(p)
end

@assert length(h5(keys)) == nMasks
close(h5)
##
h5r = h5open(h5path, "r")
h5r

# close(h5)
# imshow(h5[mask_names[1]])

## test example
h5path = joinpath("/tmp", "test.h5")
h5 = Lensman.HDF5ThreadSafe.hdf5server(h5path, "w")
@threads for i in 1:10
    h5["test$i"] = collect(i:i+5)
end
close(h5)

h5r = h5open(h5path, "r")
for k in keys(h5r)
    println("$k: ", read(h5r[k]))
end

##
h5path = joinpath("/tmp", "masks.h5")
h5r = h5open(h5path, "r")
##
for k in keys(h5r)
    println(k)
end