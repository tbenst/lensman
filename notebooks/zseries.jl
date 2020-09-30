##
import Pkg
Pkg.activate("O:\\users/tyler/code/lensman/")
using PyCall
##
prairie_view = pyimport("prairie_view")
lib = pyimport("prairie_view.lib")
##
# if error, look for popup
pv = lib.PrairieLink(env=:None)
pv.print_settings()
##
# 300 frames is wayy too slow...
# perhaps don't waste so much time at fluorphore (save ~30)
# and scan every 2um (2x speedup)
pv.max_move = 301
pv.z = -149
nFrames = 300
step_size = 1
pv.imaging_pockels = 190
pv.green_pmt = 700
pv.send("-SetZSeriesStart")
sleep(0.1)
pv.z += (nFrames-1) * step_size
pv.zseries_filename = "wholebrain"
pv.imaging_pockels = 220
pv.send("-SetZSeriesStop")
sleep(0.1)
pv.send("-SetZSeriesStepSize $step_size")
sleep(0.1)
pv.frame_averaging = 32
pv.z = 0
sleep(2)
pv.imaging_pockels = 170 # don't leave pockels high by accident


##
# TODO: python stream_frames crashes ...due to frame averaging?!?
zseries = pv.stream_frames("Zseries", nFrames, timeout=10)

##
using FileIO, RegisterQD, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob
import Base.Threads.@threads
using Unitful: μm, s

zbrain_dir = "O:\\users\\tyler\\zbrain\\"
zbrain_units = (0.798μm, 0.798μm, 2μm)

## Read Zseries from file
# tifdir = joinpath(pv.directory,"Zseries_024034.17_-000")
tifdir = raw"E:\tyler\2020-09-29_elavl3-chrmine-kv2.1_h2-6s\Zseries_024034.17_-000"
tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") & occursin("Ch3", x),
    readdir(tifdir)))
tif0 = ImageMagick.load(tiff_files[1])
H, W = size(tif0)

zseries = zeros(Normed{UInt16,16}, H, W, size(tiff_files, 1))
@threads for z in 1:size(tiff_files,1)
    zseries[:,:,z] = ImageMagick.load(tiff_files[z])
end
zseries = reverse(zseries,dims=2); # face left
size(zseries)
##
imshow(zseries)
## Read microns Per pixel (from file...)
open(glob(tifdir*"/*.xml")[1], "r") do io
    zseries_xml = read(io, String)
    global zseries_xml = xp_parse(zseries_xml)
end;
##
"apply each function to args"
funprod(functions...) = (args...)->map(x->x(args...), functions)
zeroToOne(x) = (x.-minimum(x))/(minimum(x) + maximum(x))
##
micronsPerPixel_xml = zseries_xml[xpath"""//PVStateValue[@key="micronsPerPixel"]"""][1]
# parse xml
lookup_μm = axis -> etree -> parse(Float64,
    micronsPerPixel_xml[xpath"""IndexedValue[@index="$axis"]"""][1].attr["value"])μm
microscope_units = funprod(map(lookup_μm, ["YAxis", "XAxis", "ZAxis"])...)(micronsPerPixel_xml)
# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), microscope_units)
microscope_units
## Read Zbrain H2B
h2b_zbrain = AxisArray(permutedims(
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), (0.798μm, 0.798μm, 2μm))
ps = pixelspacing(h2b_zbrain); SD = ps./maximum(ps)
size(h2b_zbrain)
## TODO: this is giving NaN...perhaps because of `centered`?
# => can copy from .ipynb that is working in notebooks...
zseries_registered = warp(zseries[:,:,196], Translation((0.,0.)),
	Base.axes(h2b_zbrain[:,:,100]));
im = RGB.(zeroToOne(h2b_zbrain[:,:,100]))
channelview(im)[2:3,:,:] .= 0
channelview(im)[2,:,:] .= zseries_registered*300
im

##
pv.close()

## 
2