# MUST RUN ON LINUX FOR ANTS support (unless compiled for windows...)

# This notebook takes a 820nm Zseries, and:
# - registers to Zbrain
# - finds best plane for imaging raphe + hab
# - finds centers of neurons
# - saves this data in case session crashes
# - creates .mat file for use with Sean's multiSLM software stims

## run cell to init VSCode Julia
true
##
import Pkg
Pkg.activate("code/lensman")
##
ENV["DISPLAY"] = "localhost:12"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots
import Base.Threads.@threads
using Unitful: μm, s
sio = pyimport("scipy.io")

zbrain_dir = ""
tifdir = ""
if Sys.iswindows()
    @error "this script only supports linux!"
else
    tifdir = "/mnt/deissero/users/tyler/b115/2020-10-12_elavl3-chrmine-kv2.1_h2b-6s_11dpf/zstack-1x/ZSeries-09172019-1403-031"
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
## Read microns Per pixel (from file...)
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

# convert to 25x units (16x is actually 14.4x, assumes Prairie View has 16x selected)
@warn "assuming 16x units but using 25x objective"
microscope_units = ((microscope_units[1:2] .* (14.4/25))..., microscope_units[3])
# microscope_units = (0.5μm, 0.5μm, 1.5μm)

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


######
# Register Zbrain to Zseries using ANTs
######

# ONLY AFFINE for speed...
moving = h2b_zbrain
fixed = zseries
tmppath = ANTsRegistration.userpath()
fixedname = joinpath(tmppath, ANTsRegistration.write_nrrd(fixed))
movingname = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
outprefix = ANTsRegistration.randstring(10)
outname = joinpath(tmppath, outprefix)
println("creating file: $outname")
# no SyN so faster, but kinda meh registration...
maxiter = 100
threshold = "1e-6"
cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
println(read(cmd, String))

affine_transform_path = glob(outprefix*"*GenericAffine.mat", tmppath)[1]
affine_transform = sio.loadmat(affine_transform_path)
# warpedname = joinpath(outname*".nrrd")
warpedname = joinpath(outname*"_Warped.nii.gz")
zbrain_registered = niread(warpedname);
println("finished registering!")
adj_zbrain_registered = adjust_histogram(imadjustintensity(zbrain_registered), Equalization())
adj_zseries = adjust_histogram(imadjustintensity(zseries), Equalization())

# visualize results as sanity check
im = RGB.(adj_zbrain_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
imshow(im)

## Read mask & transform to Zseries
# TODO: define our own hemisphere masks (split zbrain down the middle..?)
masks = matread("$zbrain_dir/MaskDatabase.mat")
raphe_mask = getMaskNameIncludes(masks, "Raphe")
# make raphe binary again
raphe_mask = antsApplyTransforms(fixedname, raphe_mask, affine_transform_path) .> 0
hab_mask = getMaskNameIncludes(masks, "Habenula")
hab_mask = antsApplyTransforms(fixedname, hab_mask, affine_transform_path) .> 0

raphe_voxels_by_plane = mapslices(sum, raphe_mask, dims=[1,2])
max_raphe_voxels_in_plane = maximum(raphe_voxels_by_plane)
hab_voxels_by_plane = mapslices(sum, hab_mask, dims=[1,2])
max_hab_voxels_in_plane = maximum(hab_voxels_by_plane)

# find plane with large fraction of habenula + raphe
planeToImage = argmax(minimum.(
    zip( vec(raphe_voxels_by_plane ./ max_raphe_voxels_in_plane),
         vec(hab_voxels_by_plane ./ max_hab_voxels_in_plane)
        )))
@show planeToImage

## Segmentation

# CREATE NUCLEUS TEMPLATE
nucleus_r = (2.75μm/microscope_units[1])μm
space = 0.3μm # WARN: need to adjust CirclePointRadius, too
cell_r = nucleus_r + space
template_size = Int32(ceil(2*(cell_r)/microscope_units[1]))
# ensure odd size
template_size += abs((template_size) % 2 - 1)
# we make large and later resize for anti-aliasing
nucleus_template = Gray.(zeros(513,513))
draw!(nucleus_template, Ellipse(CirclePointRadius(256,256,192)))
nucleus_template = Float64.(nucleus_template)
# nucleus_template[nucleus_template.==0.] .= -1.
sz = (template_size, template_size)
σ = map((o,n)->0.75*o/n, size(nucleus_template), sz)
kern = KernelFactors.gaussian(σ)

nucleus_template = imresize(imfilter(nucleus_template, kern, NA()), sz)

# RESCALE IMAGE
# rescaled_img = adjust_gamma(imadjustintensity(img), 0.45)
img = zseries[:,:,planeToImage]
rescaled_img = (adjust_histogram(imadjustintensity(img), Equalization()))
# mask out background
laplace_mask = opening(imadjustintensity(morpholaplace(rescaled_img)).<0.5)
Gray.(laplace_mask)
# img = rescaled_img .* laplace_maskf1.
