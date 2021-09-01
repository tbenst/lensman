# unfinished, but has improved segmentation I think..?
##
true
##
import Pkg
Pkg.activate(".")
@warn raw"ensure using 25x AssignRoot in Z:\mSLM_B115\Main\multiSLM_Server\Startup, else SLM will be off"
@warn raw"ensure using Olympus 25x in Prairie View, else ETL will be off"
ENV["DISPLAY"] = "localhost:11"
##
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates
    import Base.Threads.@threads
using Unitful: μm, m, s
sio = pyimport("scipy.io")
morph = pyimport("skimage.morphology")

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

zbrain_dir = ""
tifdir = ""
if Sys.iswindows()
    @error "this script only supports linux!"
else
    # first we 
    tifdir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish3/ZSeries-structural-052"
    zbrain_dir = "/mnt/deissero/users/tyler/zbrain"
end
fishDir, _ = splitdir(tifdir)

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
size(zseries)

#####
    ## Read microns Per pixel (from file...)
#####

# open(glob("*.xml",tifdir)[1], "r") do io
#     zseries_xml = read(io, String)
#     global zseries_xml = xp_parse(zseries_xml)
# end

# "apply each function to args"
# funprod(functions...) = (args...)->map(x->x(args...), functions)

# micronsPerPixel_xml = zseries_xml[xpath"""//PVStateValue[@key="micronsPerPixel"]"""][1]
# # parse xml
# lookup_μm = axis -> etree -> parse(Float64,
#     micronsPerPixel_xml[xpath"""IndexedValue[@index="$axis"]"""][1].attr["value"])μm
# microscope_units = funprod(map(lookup_μm, ["YAxis", "XAxis", "ZAxis"])...)(micronsPerPixel_xml)

# convert to 25x units (16x is actually 14.4x, assumes Prairie View has 16x selected)
# @info "assuming 16x units but using 25x objective"
# microscope_units = ((microscope_units[1:2] .* (14.4/25))..., microscope_units[3])
# microscope_units = (0.5μm, 0.5μm, 1.5μm)

# TODO: make Olympus 25x in bruker have proper units....
@info "hardcoded microscope units"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)

# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), microscope_units)
adj_zseries = adjust_histogram(imadjustintensity(zseries), GammaCorrection(0.5))
adj_zseries = AxisArray(adj_zseries, (:y, :x, :z), microscope_units)
@assert microscope_units[1] > 0.1μm
microscope_units

#####
# Read Zbrain H2B
#####
function flipZbrain(x)
    x = reverse(x,dims=2); # face right
    x = reverse(x,dims=3); # top to bottom
    x
end


h2b_zbrain = AxisArray(permutedims(
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), zbrain_units)
h2b_zbrain = reverse(h2b_zbrain,dims=2); # face right
h2b_zbrain = reverse(h2b_zbrain,dims=3); # top to bottom
ps = pixelspacing(h2b_zbrain); SD = ps./maximum(ps)
adj_h2b_zbrain = adjust_histogram(imadjustintensity(h2b_zbrain), Equalization())
adj_h2b_zbrain = AxisArray(adj_h2b_zbrain, (:y, :x, :z), zbrain_units)
size(h2b_zbrain)
# ## Examine our registration starting point...
# im = RGB.(adjust_histogram(imadjustintensity(h2b_zbrain), Equalization()))
# channelview(im)[2:3,:,:,:] .= 0
# # channelview(im)[2,:,:,:] .= adj_zseries
# imshow(im)
# # warp(zseries[:,:,196], Translation((0.,0.))
#


######
## Register Zbrain to Zseries using ANTs
######

# ONLY AFFINE for speed...
# we use histogram matching in ANTs (although maybe our method is better...?)
fixed = zseries
moving = h2b_zbrain
fixedname = joinpath(fishDir, "zseries.nrrd")
save(fixedname, fixed)
movingname = joinpath(fishDir, "zbrain_h2b.nrrd")
save(movingname, moving)
outprefix = Dates.format(DateTime(now()), DateFormat("YmmddTHHMMSSsZ"))
outname = joinpath(fishDir, outprefix)
println("creating files with name: $outname")
interpolation = "BSpline"
histmatch = 1
maxiter = 200
threshold = "1e-6"
initial_moving_type = 1
sampling_frac = 1.0
synThreshold = "1e-6"
# takes ~5 min
cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation $interpolation --use-histogram-matching $histmatch -r \[$fixedname, $movingname,$initial_moving_type\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,    $sampling_frac\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,$sampling_frac\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
println("calling ANTs...")
@time println(read(cmd, String))

# affine_transform_path = glob(outprefix*"_SyN"*"*GenericAffine.mat", fishDir)[1]
affine_transform_path = glob(outprefix*"*GenericAffine.mat", fishDir)[1]
affine_transform = sio.loadmat(affine_transform_path)
# warpedname = joinpath(outname*".nrrd")
warpedname = joinpath(outname*"_Warped.nii.gz")

# if this fails, tmpfs may be out of space
# try `rm -rf /run/user/1000/tyler/ANTs/`
zbrain_registered = niread(warpedname);

println("finished registering!")
adj_zbrain_registered = adjust_histogram(imadjustintensity(zbrain_registered), Equalization());

## visualize results as sanity check
im = RGB.(adj_zbrain_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
imshow(im)

## Read mask & transform to Zseries
# TODO: define our own hemisphere masks (split zbrain down the middle..?)
masks = matread("$zbrain_dir/MaskDatabase.mat")
raphe_mask = flipZbrain(getMaskNameIncludes(masks, "Raphe"))
# make raphe binary again
raphe_mask = antsApplyTransforms(fixedname, raphe_mask, affine_transform_path) .> 0
hab_mask = flipZbrain(getMaskNameIncludes(masks, "Habenula"))
hab_mask = antsApplyTransforms(fixedname, hab_mask, affine_transform_path) .> 0

eyes_mask = flipZbrain(getMaskNameIncludes(masks, "Eyes"))
eyes_mask = antsApplyTransforms(fixedname, eyes_mask, affine_transform_path) .> 0

## neuron locations by plane, takes ~30 seconds
neuronCenterMask = zeros(Bool, size(zseries)...)
for z in 1:size(neuronCenterMask, 3)
    # neuronCenterMask[:,:,z] = findNeurons(zseries[:,:,z])
    # quick debugging hack as Lensman not compiling..?
    neuronCenterMask[:,:,z] = fN(zseries[:,:,z])
end

raphe_nNeurons_by_plane = mapslices(sum, raphe_mask .& neuronCenterMask, dims=[1,2])
max_raphe_nNeurons_in_plane = maximum(raphe_nNeurons_by_plane)
hab_nNeurons_by_plane = mapslices(sum, hab_mask .& neuronCenterMask, dims=[1,2])
max_hab_nNeurons_in_plane = maximum(hab_nNeurons_by_plane)
# auto determine plane
planeHab = Tuple(argmax(hab_nNeurons_by_plane))[3]
planeRaphe = Tuple(argmax(raphe_nNeurons_by_plane))[3]

step = Int(round((planeRaphe - planeHab)//4))

# we choose only planes from zseries for easy segmentation
planesToImage = planeHab .+ collect(0:4) .* step
zeroPlane = planesToImage[3]
zseries_step = 2. # microns
##
nNeuronPerPlane = 25
targetLocs = [sample(findall((~).(eyes_mask[:,:,planesToImage[z]]) .& neuronCenterMask[:,:,planesToImage[z]]),
    nNeuronPerPlane, replace=false) for z in 1:5]
arrCart2arr(x) = Float64.(hcat(collect.(Tuple.(x))...))
targetsByPlane =  arrCart2arr.(targetLocs)
targetsZoffset = zseries_step .* (zeroPlane .- planesToImage) .* 10^-6
targetsByPlane = [vcat(targetsByPlane[z],
    ones(Float64,1,nNeuronPerPlane)*targetsZoffset[z])'
    for z in 1:5]

for z in 1:5
    # account for sean's code using 512 x 512 space
    # while we image in 1024x1024
    targetsByPlane[z][:,[1,2]] ./= 2
    # swap for sean's ordering convention
    targetsByPlane[z][:,[2,1,3]]
end

## Visualize stim targets
im = RGB.(adj_zseries[:,:,planesToImage])
stim_points = zeros(Bool,size(im))
for z in 1:5
    stim_points[targetLocs[z],z] .= true
    stim_points[:,:,z] = dilate(dilate(stim_points[:,:,z]))
end
channelview(im)[[1,3],:,:,:] .= 0
channelview(im)[1,:,:,:] .= 0.5*float(stim_points)
imshow(im)


##
im = adjust_histogram(imadjustintensity(h2b_zbrain), GammaCorrection(0.5))
imshow(im)
##
imshow(h2b_zbrain)
##
im = adjust_histogram(imadjustintensity(h2b_zbrain), Equalization())
imshow(im)


##
plane = zseries[:,:,planeHab]
plane = adjust_histogram(imadjustintensity(plane), GammaCorrection(0.5))
threshold = otsu_threshold(plane)
foreground = plane .> threshold
feat_img = match_template(Float64.(foreground), disk(4), pad_input=true, mode="constant", constant_values=0)
local_maxi = peak_local_max(feat_img, indices=false, footprint=disk(4), exclude_border=false, threshold_rel=.1)
# Gray.(dilate(dilate(local_maxi)))
seeds = label_components(local_maxi)
segments = morph.watershed(-feat_img, seeds, mask=foreground)
# segments = watershed(Float64.(-feat_img), seeds; compactness=1e1, mask=foreground)

# colored_segments = map(i->get_random_color(i), labels_map(segments))
colored_segments = map(i->get_random_color(i), segments)
imshow(colored_segments)

##
min = 5
max = 9*9
valid_segments = [k for (k,v) in segment_pixel_count(segments)
    if (v>=min) & (v<=max)]

im = RGB.(plane)
channelview(im)[[1,3],:,:] .= 0
channelview(im)[1,:,:]  .= in.(labels_map(segments), Ref(Set(valid_segments)))
im
##
# laplace_mask = opening(imadjustintensity(morpholaplace(plane)).<0.5)
# Gray.(laplace_mask)
# segment_nuclei(channelview(plane))

##
# TODO: switch to GLMakie
# but need to add libglfw3 to nix shell
using Makie, GeometryTypes, Colors, GLMakie, AbstractPlotting

img = rand(100,100)
scene = Scene()
points = [Makie.Point2f0(cos(t), sin(t)) for t in LinRange(0, 2pi, 20)]
colors = 1:20
Makie.scatter!(scene, points, color = colors, markersize = 15)]