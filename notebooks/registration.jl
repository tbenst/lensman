# file to load zseries from Oak, register to z-brain using ANTs
# and finally segment cells using a watershed

ENV["DISPLAY"] = "localhost:11.0"
using FileIO, RegisterQD, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation
import Base.Threads.@threads
using Unitful: μm, s

zbrain_dir = ""
tifdir = ""
if Sys.iswindows()
    tifdir = raw"E:\tyler\2020-10-05_elavl3-chrmine-kv2.1_h2b-6s\ZSeries-25X-default-exponential-zpower-64avg-034"
    zbrain_dir = "O:\\users\\tyler\\zbrain\\"
else
    tifdir = "/mnt/deissero/users/tyler/b115/2020-09-29_elavl3-chrmine-kv2.1_h2-6s/Zseries_024034.17_-000"
    zbrain_dir = "/mnt/deissero/users/tyler/zbrain"
end
zbrain_units = (0.798μm, 0.798μm, 2μm)

# Read Zseries from file
# tifdir = joinpath(pv.directory,"Zseries_024034.17_-000")
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
#####
# Read microns Per pixel (from file...)
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

@warn "temp hack for bad units"
microscope_units = (0.5μm, 0.5μm, 1.5μm)

# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), microscope_units)
microscope_units
#####
# Read Zbrain H2B
#####
h2b_zbrain = AxisArray(permutedims(
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), (0.798μm, 0.798μm, 2μm))
ps = pixelspacing(h2b_zbrain); SD = ps./maximum(ps)
size(h2b_zbrain)

#####
# Inital transform guess
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
# init_x_translation = -50μm/x_units
init_x_translation = -70μm/x_units
init_z_translation = -40

# 3D preview pre-affine
initial_tfm = Translation(init_y_translation, init_x_translation, init_z_translation)
initial_affine = AffineMap([1/y_scaling 0. 0.;
                            0. 1/x_scaling 0.;
                            0. 0. 1/z_scaling], initial_tfm.translation)
replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
zseries_registered = replace_nan.(warp(zseries, initial_affine,
    Base.axes(h2b_zbrain)));
size(zseries_registered)
#
im = RGB.(zeroToOne(h2b_zbrain))
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= zseries_registered*300
imshow(im)

##########
## ANTs
##########
stageaff = Stage(h2b_zbrain, Global("Affine"))
stagesyn = Stage(h2b_zbrain, SyN())
imgw = register(h2b_zbrain, zseries_registered, [stageaff, stagesyn])
# crashes on reading in nrrd due to zlib.so not found :/ 
# my LD_LIBRARY_PATH is effed..... due to conda shell
# replicate crash
# regular juila-env works!
imshow(imgw)
##
# show registration onto our brain...
im = RGB.(zeroToOne(h2b_zbrain))
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= imgw*300
imshow(im)

###############
## SEGMENTATION
###############
##
threshold = otsu_threshold(float.(h2b_zbrain))
##
imshow(h2b_zbrain.>threshold)
##
seg = fast_scanning(h2b_zbrain[:,:,100],threshold)
seg
##
using TestImages
img = testimage("cameraman")
threshold = convert(Float64, otsu_threshold(img))
seg = fast_scanning(img,0.1)
imshow(labels_map(seg))
##########
## ARCHIVE
##########

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
imshow(zseries)
##

## BASELINE
zseries_registered = warp(zseries[:,:,196], Translation((0.,0.)),
	Base.axes(h2b_zbrain[:,:,100]));
im = RGB.(zeroToOne(h2b_zbrain[:,:,100]))
channelview(im)[2:3,:,:] .= 0
channelview(im)[2,:,:] .= zseries_registered*300
im

##