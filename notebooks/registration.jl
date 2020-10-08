# file to load zseries from Oak, register to z-brain using ANTs
# and finally segment cells using a watershed
# how can we use PyCall but not call this binary...?
# why aren't we in a conda-shell :/
##
true
##
ENV["DISPLAY"] = "localhost:11.0"
using FileIO, RegisterQD, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI
import Base.Threads.@threads
using Unitful: μm, s
sio = pyimport("scipy.io")

zbrain_dir = ""
tifdir = ""
if Sys.iswindows()
    tifdir = raw"E:\tyler\2020-10-05_elavl3-chrmine-kv2.1_h2b-6s\ZSeries-25X-default-exponential-zpower-64avg-034"
    zbrain_dir = "O:\\users\\tyler\\zbrain\\"
else
    tifdir = "/mnt/deissero/users/tyler/b115/2020-10-05_elavl3-chrmine-kv2.1_h2b-6s/ZSeries-25X-default-exponential-zpower-64avg-034"
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
im = RGB.(imadjustintensity(h2b_zbrain))
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= zseries_registered*300
# imshow(im)
true
##########
## ANTs
# parameters from https://academic.oup.com/gigascience/article/6/8/gix056/4085311#supplementary-data
function marquart_register(fixed::AbstractArray, moving::AbstractArray)
    fixedname = ANTsRegistration.write_nrrd(fixed)
    movingname = ANTsRegistration.write_nrrd(moving)
    outname = ANTsRegistration.randstring(10)
    cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t SyN\[0.05,6,0.5\] -m CC\[$fixedname, $movingname,1,2\] -c \[200x200x200x200x10,1e-7,10\] --shrink-factors 12x8x4x2x1 --smoothing-sigmas 4x3x2x1x0vox`
    # cmd = `antsRegistration -h`
    read(cmd, String)
end

fixed = h2b_zbrain
moving = zseries
tmppath = ANTsRegistration.userpath()
fixedname = joinpath(tmppath, ANTsRegistration.write_nrrd(fixed))
movingname = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
outprefix = ANTsRegistration.randstring(10)
outname = joinpath(tmppath, outprefix)
# no SyN so faster, but kinda meh registration...
cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
    # cmd = `antsRegistration -h`
println(read(cmd, String))

# `antsRegistration -d 3 –float 1 -o [fish1_, fish1_Warped.nii.gz] –interpolation WelchWindowedSinc –use-histogram-matching 0 -r [ref/vglut-ref.nii, fish1–01.nii.gz,1] -t rigid[0.1] -m MI[ref/vglut-ref.nii, fish1–01.nii.gz,1,32, Regular,0.25] -c [200 × 200 × 200 × 0,1e-8,10] –shrink-factors 12 × 8 × 4 × 2 –smoothing-sigmas 4 × 3 × 2 × 1vox -t Affine[0.1] -m MI[ref/vglut-ref.nii, fish1–01.nii.gz,1,32, Regular,0.25] -c [200 × 200 × 200 × 0,1e-8,10] –shrink-factors 12 × 8 × 4 × 2 –smoothing-sigmas 4 × 3 × 2 × 1vox -t SyN[0.05,6,0.5] -m CC[ref/vglut-ref.nii, fish1–01.nii.gz,1,2] -c [200 × 200 × 200 × 200 × 10,1e-7,10] –shrink-factors 12 × 8 × 4 × 2 × 1 –smoothing-sigmas 4 × 3 × 2 × 1 × 0vox `

# ERROR:  The number of smoothing sigmas does not match the number of levels.
# println(marquart_register(h2b_zbrain, zseries_registered))
# println(fast_register(h2b_zbrain, zseries_registered))
##
affine_transform_path = glob(outprefix*"*GenericAffine.mat", tmppath)[1]
affine_transform = sio.loadmat(affine_transform_path)
# warpedname = joinpath(outname*".nrrd")
warpedname = joinpath(outname*"_Warped.nii.gz")
# TODO: make NIfTI use CodecZlib instead...
registered = niread(warpedname)

im = RGB.(imadjustintensity(h2b_zbrain))
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= imadjustintensity(registered)
imshow(im)

##
2*4
##



stageaff = Stage(h2b_zbrain, Global("Affine"))
stagesyn = Stage(h2b_zbrain, SyN())
up = ANTsRegistration.userpath()
outname = ANTsRegistration.randstring(10)
tfmname = joinpath(up, outname*"_warp")
warpedname = joinpath(up, outname*".nrrd")
output = (tfmname, warpedname)


register(output, h2b_zbrain, zseries_registered,
    [stageaff, stagesyn])
    # [stageaff])
##
# TODO: is this v4 .mat file...?
# affine_transform = matread(glob(outname*"*GenericAffine.mat", up)[1])
affine_transform = sio.loadmat(glob(outname*"*GenericAffine.mat", up)[1])
##
run(`antsApplyTransforms -h`)
# TODO: use antsApplyTransforms to tranfrom masks to our fish
##
imgw = load(warpedname)
# TODO: warp 
# crashes on reading in nrrd due to zlib.so not found :/ 
# my LD_LIBRARY_PATH is effed..... due to conda shell
# replicate crash
# regular juila-env works!
# imshow(imgw)
##
# show registration onto our brain...
im = RGB.(imadjustintensity(h2b_zbrain))
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
zseries_affine = qd_affine(imadjustintensity(h2b_zbrain),
	imadjustintensity(zseries), [100,100, 50];
    initial_tfm=initial_affine)
##
zseries_registered = warp(zseries, zseries_affine[1],
	Base.axes(h2b_zbrain));
im = RGB.(imadjustintensity(h2b_zbrain))
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
im = RGB.(imadjustintensity(h2b_zbrain[:,:,100]))
channelview(im)[2:3,:,:] .= 0
channelview(im)[2,:,:] .= zseries_registered*300
im

##