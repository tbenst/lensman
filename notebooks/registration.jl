# file to load zseries from Oak, register to z-brain using ANTs
# and finally segment cells using a watershed
# how can we use PyCall but not call this binary...?
# why aren't we in a conda-shell :/
##
true
##
# ENV["DISPLAY"] = "localhost:11.0" # remote-ssh
ENV["DISPLAY"] = ":0"
using FileIO, RegisterQD, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment
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

##########
# ANTs
# parameters from https://academic.oup.com/gigascience/article/6/8/gix056/4085311#supplementary-data
function marquart_register(fixed::AbstractArray, moving::AbstractArray)
    fixedname = ANTsRegistration.write_nrrd(fixed)
    movingname = ANTsRegistration.write_nrrd(moving)
    outname = ANTsRegistration.randstring(10)
    cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t SyN\[0.05,6,0.5\] -m CC\[$fixedname, $movingname,1,2\] -c \[200x200x200x200x10,1e-7,10\] --shrink-factors 12x8x4x2x1 --smoothing-sigmas 4x3x2x1x0vox`
    # cmd = `antsRegistration -h`
    read(cmd, String)
end


# ONLY AFFINE for speed...
# fixed = h2b_zbrain
# moving = zseries
moving = h2b_zbrain
fixed = zseries
tmppath = ANTsRegistration.userpath()
fixedname = joinpath(tmppath, ANTsRegistration.write_nrrd(fixed))
movingname = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
outprefix = ANTsRegistration.randstring(10)
outname = joinpath(tmppath, outprefix)
println("creating file: $outname")
# no SyN so faster, but kinda meh registration...
# cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[200x200x200x0,1e-8,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
maxiter = 200
maxiter = 100
threshold = "1e-8"
threshold = "1e-6"
cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
# cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
    # cmd = `antsRegistration -h`
println(read(cmd, String))

affine_transform_path = glob(outprefix*"*GenericAffine.mat", tmppath)[1]
affine_transform = sio.loadmat(affine_transform_path)
# warpedname = joinpath(outname*".nrrd")
warpedname = joinpath(outname*"_Warped.nii.gz")
registered = niread(warpedname);
println("finished registering!")
adj_registered = adjust_histogram(imadjustintensity(registered), Equalization())
adj_zseries = adjust_histogram(imadjustintensity(zseries), Equalization())

im = RGB.(adj_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
imshow(im)
##
# this works!! actually labels raphe decently well
zseries_in_zbrain_space = niread("/tmp/tyler/ANTs/eHDo5NtK12_WarpedInv.nii.gz")
im = RGB.(adjust_histogram(imadjustintensity(zseries_in_zbrain_space), Equalization()))
channelview(im)[1,:,:,:] .= 0
channelview(im)[3,:,:,:] .= Gray.(float32.(raphe_mask))
imshow(im)

##
########
## WIP
########
masks = matread("$zbrain_dir/MaskDatabase.mat")
##
"Read (sparse) mask from Zbrain atlas and reshape to dense Array."
function read_mask(masks, idx; W=621, H=1406, Z=138, units=zbrain_units)
    example_mask = reshape(Array(masks["MaskDatabase"][:,idx]), (H, W, Z))
    # reshape as [621, 1406, 138]
    example_mask = permutedims(example_mask, (2,1,3))
    AxisArray(example_mask, (:y, :x, :z), units)
end

##
masks_raphe_idxs = findall(occursin.("Raphe", masks["MaskDatabaseNames"])) |>
    x -> map(ci -> ci[2], x)
raphe_mask = foldl((a,n)->a .| read_mask(masks, n), masks_raphe_idxs;
    init=read_mask(masks,masks_raphe_idxs[1]))
imshow(Gray.(raphe_mask))

## transform to our brain...
maskmovingname = joinpath(tmppath, ANTsRegistration.write_nrrd(
    AxisArray(float32.(raphe_mask), (:y, :x, :z), zbrain_units)))
# maskmovingname = joinpath(tmppath, ANTsRegistration.write_nrrd(h2b_zbrain))
# maskoutname = joinpath(tmppath, ANTsRegistration.randstring(10)*".nii.gz")
run(`antsApplyTransforms --float -d 3 -i $maskmovingname -r $fixedname -t $affine_transform_path -o $maskoutname`)
# run(`antsApplyTransforms --float -d 3 -i $fixedname -r $movingname -t \[ $affine_transform_path, 1 \] -o $maskoutname`)
# run(`antsApplyTransforms --float -d 3 -i $maskmovingname -r $fixedname -n NearestNeighbor -t $warpedname -o $maskoutname`)
# run(`antsApplyTransforms --float -d 3 -i $maskmovingname -r $fixedname -n NearestNeighbor -t $affine_transform_path -t $warpedname -o $maskoutname`)

zseries_raphe_mask = niread(maskoutname)
sum(zseries_raphe_mask)

imshow(zseries_raphe_mask)
##
imshow(adj_registered)
##
im = RGB.(adj_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
channelview(im)[3,:,:,:] .= zseries_raphe_mask
imshow(im)
##
h2b_zbrain_adj = adjust_histogram(imadjustintensity(h2b_zbrain), Equalization())
im = RGB.(h2b_zbrain_adj)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[3,:,:,:] .= Gray.(float32.(raphe_mask))
imshow(im)
#########
## TODO
#########


## compare to qd_affine
zseries_affine = qd_affine(zeroToOne(h2b_zbrain[80,:,:]),
	zeroToOne(zseries[196,:,:]), [50,50];
    initial_tfm=initial_affine)


##
2*4
##
##


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
