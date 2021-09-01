##
# import Pkg
# Pkg.activate("O:\\users/tyler/code/lensman/")
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
ENV["DISPLAY"] = "localhost:11.0"
using FileIO, RegisterQD, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration
import Base.Threads.@threads
using Unitful: μm, s

# zbrain_dir = "O:\\users\\tyler\\zbrain\\"
zbrain_dir = "/mnt/deissero/users/tyler/zbrain"
zbrain_units = (0.798μm, 0.798μm, 2μm)

## Read Zseries from file
# tifdir = joinpath(pv.directory,"Zseries_024034.17_-000")
# tifdir = raw"E:\tyler\2020-09-29_elavl3-chrmine-kv2.1_h2-6s\Zseries_024034.17_-000"
tifdir = "/mnt/deissero/users/tyler/b115/2020-09-29_elavl3-chrmine-kv2.1_h2-6s/Zseries_024034.17_-000"
tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") & occursin("Ch3", x),
    readdir(tifdir)))
tif0 = ImageMagick.load(tiff_files[1])
H, W = size(tif0)

zseries = zeros(Normed{UInt16,16}, H, W, size(tiff_files, 1))
@threads for z in 1:size(tiff_files,1)
    zseries[:,:,z] = ImageMagick.load(tiff_files[z])
end
zseries = reverse(zseries,dims=2); # face left
zseries = reverse(zseries,dims=3); # bottom to top
size(zseries)
##
imshow(zseries)
#####
## Read microns Per pixel (from file...)
#####
open(glob("*.xml",tifdir)[1], "r") do io
    zseries_xml = read(io, String)
    global zseries_xml = xp_parse(zseries_xml)
end;

"apply each function to args"
funprod(functions...) = (args...)->map(x->x(args...), functions)
zeroToOne(x) = (x.-minimum(x))/(minimum(x) + maximum(x))

micronsPerPixel_xml = zseries_xml[xpath"""//PVStateValue[@key="micronsPerPixel"]"""][1]
# parse xml
lookup_μm = axis -> etree -> parse(Float64,
    micronsPerPixel_xml[xpath"""IndexedValue[@index="$axis"]"""][1].attr["value"])μm
microscope_units = funprod(map(lookup_μm, ["YAxis", "XAxis", "ZAxis"])...)(micronsPerPixel_xml)
# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), microscope_units)
microscope_units
#####
## Read Zbrain H2B
#####
h2b_zbrain = AxisArray(permutedims(
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), (0.798μm, 0.798μm, 2μm))
ps = pixelspacing(h2b_zbrain); SD = ps./maximum(ps)
size(h2b_zbrain)


## BASELINE
zseries_registered = warp(zseries[:,:,196], Translation((0.,0.)),
	Base.axes(h2b_zbrain[:,:,100]));
im = RGB.(zeroToOne(h2b_zbrain[:,:,100]))
channelview(im)[2:3,:,:] .= 0
channelview(im)[2,:,:] .= zseries_registered*300
im

#####
## Inital transform guess
#####
zbrain_x_units = h2b_zbrain[Axis{:x}][2] - h2b_zbrain[Axis{:x}][1]
zbrain_y_units = h2b_zbrain[Axis{:y}][2] - h2b_zbrain[Axis{:y}][1]
zbrain_z_units = h2b_zbrain[Axis{:z}][2] - h2b_zbrain[Axis{:z}][1]
x_units = zseries[Axis{:x}][2] - zseries[Axis{:x}][1]
y_units = zseries[Axis{:y}][2] - zseries[Axis{:y}][1]
z_units = zseries[Axis{:z}][2] - zseries[Axis{:z}][1]
x_scaling = x_units/zbrain_x_units
y_scaling = y_units/zbrain_y_units
z_scaling = z_units/zbrain_z_units
(zbrain_x_units, zbrain_y_units, zbrain_z_units), (z_units, y_units, x_units)

init_y_translation = -(size(h2b_zbrain[Axis{:y}])[1] - size(zseries[Axis{:y}])[1]*y_scaling)/2
ignore_zbrain_x = 
# init_x_translation = -(size(h2b_zbrain[Axis{:x}])[1] - size(zseries[Axis{:x}])[1]*x_scaling)/2 # TODO: thi
init_x_translation = -50μm/x_units
initial_tfm = Translation(init_y_translation, init_x_translation)
initial_affine = AffineMap([1/y_scaling 0.; 0. 1/x_scaling], initial_tfm.translation)

# initial_affine = AffineMap([2 0.; 0. 2], [-100, -500])

zseries_registered = warp(zseries[:,:,196], initial_affine,
	Base.axes(h2b_zbrain[:,:,80]));
im = RGB.(zeroToOne(h2b_zbrain[:,:,67]))
channelview(im)[2:3,:,:] .= 0
channelview(im)[2,:,:] .= zseries_registered*300
im

## AFFINE REGISTER
zseries_affine = qd_affine(zeroToOne(h2b_zbrain[80,:,:]),
	zeroToOne(zseries[196,:,:]), [50,50];
    initial_tfm=initial_affine)

## 3D preview pre-affine
initial_tfm = Translation(init_y_translation, init_x_translation, 0)
initial_affine = AffineMap([1/y_scaling 0. 0.;
                            0. 1/x_scaling 0.;
                            0. 0. 1/z_scaling], initial_tfm.translation)
zseries_registered = warp(zseries, initial_affine,
	Base.axes(h2b_zbrain));
im = RGB.(zeroToOne(h2b_zbrain))
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= zseries_registered*300
imshow(im)
## AFFINE REGISTER
zseries_affine = qd_affine(zeroToOne(h2b_zbrain),
	zeroToOne(zseries), [100,100, 50];
    initial_tfm=initial_affine)
##
zseries_registered = warp(zseries, zseries_affine[1],
	Base.axes(h2b_zbrain));
im = RGB.(zeroToOne(h2b_zbrain))
channelview(im)[2 .= 0
channelview(im) .= zseries_registered*300
im
##


##
pv.close()

## 
2