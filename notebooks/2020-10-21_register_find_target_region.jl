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
Pkg.activate(".")
##
ENV["DISPLAY"] = "localhost:11"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots
    import Base.Threads.@threads
using Unitful: μm, s
sio = pyimport("scipy.io")

zbrain_dir = ""
tifdir = ""
if Sys.iswindows()
    @error "this script only supports linux!"
else
    # tifdir = "/mnt/deissero/users/tyler/b115/2020-10-05_elavl3-chrmine-kv2.1_h2b-6s/ZSeries-25X-default-exponential-zpower-64avg-034/"
    # tifdir = "/mnt/deissero/users/tyler/b115/2020-10-20_elavl3-chrmine-kv2.1_h2b6s_6dpf/ZSeries-12112019-1026-002"
    # tifdir = "/mnt/deissero/users/tyler/b115/2020-10-21_elavl3-chrmine-kv2.1_h2b6s_dpf/ZSeries-12112019-1026-004"
    # tifdir = "/mnt/deissero/users/tyler/b115/2020-10-21_elavl3-chrmine-kv2.1_h2b6s_dpf/fish2/ZSeries-09172019-1403-036"
    tifdir = "/mnt/deissero/users/tyler/b115/2020-10-22_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/ZSeries-09172019-1403-039"
    tifdir = "/mnt/deissero/users/tyler/b115/2020-10-22_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/ZSeries-09172019-1403-041"

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
# TODO: we should reverse atlas NOT our imaging to make creating .mat easier...
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

# convert to 25x units (16x is actually 14.4x, assumes Prairie View has 16x selected)
# @warn "assuming 16x units but using 25x objective"
# microscope_units = ((microscope_units[1:2] .* (14.4/25))..., microscope_units[3])
# microscope_units = (0.5μm, 0.5μm, 1.5μm)

# TODO: make Olympus 25x in bruker have proper units....
@warn "hardcoded microscope units"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)

# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), microscope_units)
adj_zseries = adjust_histogram(imadjustintensity(zseries), Equalization())
adj_zseries = AxisArray(adj_zseries, (:y, :x, :z), microscope_units)
@assert microscope_units[1] > 0.1μm
microscope_units

#####
# Read Zbrain H2B
#####
h2b_zbrain = AxisArray(permutedims(
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), zbrain_units)
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
# moving = adj_h2b_zbrain
# fixed = adj_zseries
tmppath = ANTsRegistration.userpath()
fixedname = joinpath(tmppath, ANTsRegistration.write_nrrd(fixed))
movingname = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
outprefix = ANTsRegistration.randstring(10)
outname = joinpath(tmppath, outprefix)
# outname = joinpath(tmppath, outprefix*"_SyN")
println("creating file: $outname")
# no SyN so faster, but kinda meh registration...
maxiter = 100
threshold = "1e-6"
synThreshold = "1e-5"
# originally was -r instead of -q
# maybe try both?
# -q is --initial-fixed-transform
# -r is --initial-moving-transform
# no SyN
cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,0\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
# SyN - takes ~20 minutes, TODO can this be sped up??
# cmd = `antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation WelchWindowedSinc --use-histogram-matching 0 -r \[$fixedname, $movingname,1\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,0.25\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t SyN\[0.05,6,0.5\] -m CC\[$fixedname, $movingname,1,2\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$synThreshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
println(read(cmd, String))

# affine_transform_path = glob(outprefix*"_SyN"*"*GenericAffine.mat", tmppath)[1]
affine_transform_path = glob(outprefix*"*GenericAffine.mat", tmppath)[1]
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
        
planeHab = Tuple(argmax(hab_voxels_by_plane))[3]
planeRaphe = Tuple(argmax(raphe_voxels_by_plane))[3]
zero_plane = Int(round((planeRaphe+planeHab)//2))

@warn "assuming each plane is 2μm apart"
zseries_step_size = 2
step_size = (planeHab - zero_plane) * zseries_step_size
numZplanes = size(zseries,3)
# note that our zseries is inverted from PV...
# TODO: flip this so more intuitive...?
# TODO: why are we flipped lr...?
# println("Set Z 'zero' to be plane $(numZplanes-zero_plane) and step size $step_size μm to image 3 planes")
# we move up one plane (2um) to account for spectral shift
println("Set first plane to  $(numZplanes-planeHab-1) and step size $step_size μm to image 3 planes")
# TODO: read state from xml file and communicate with PV so we get zseries microns correct...?

println("showing prairie view plane numbers")
@show numZplanes - zero_plane
@show numZplanes - planeHab
@show numZplanes - planeRaphe

##
im = RGB.(raphe_mask)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
channelview(im)[3,:,:,:] .= hab_mask
imshow(im)


## Segmentation
local_maxi = findNeurons(zseries[:,:,planeHab])
@show sum(local_maxi)
maxi_mask = (dilate ∘ dilate ∘ dilate)(local_maxi) .> 0
Gray.(maxi_mask)

## MAKE PROTOCOL
# we image 3 planes: best hab, best raphe, and best compromise
# we have three stim groups:
#   (1) left hab (best hab plane) 
#   (2) right hab (best hab plane) 
#   (3) raphe (best raphe plane) 
plane_neurons = findNeurons(zseries[:,:,planeHab])
right_mask = zeros(Bool, size(plane_neurons)...)
dividing_line = Int(round(mean(mask2ind(raphe_mask)[:,1])))
@warn "dividing_line set to $dividing_line, check if reasonable?"
right_mask[1:dividing_line,:] .= 1
left_mask = zeros(Bool, size(plane_neurons)...)
left_mask[dividing_line+1:end,:] .= 1
hab_neurons = (plane_neurons .& hab_mask[:,:,planeHab])
right_hab_neurons = hab_neurons .& right_mask
left_hab_neurons = hab_neurons .& left_mask
@show sum(right_hab_neurons)
@show sum(left_hab_neurons)

plane_neurons = findNeurons(zseries[:,:,planeRaphe])
raphe_neurons = plane_neurons .& raphe_mask[:,:,planeRaphe]
@show sum(raphe_neurons)

## sample consistent number of neurons
nNeuronsForStim = 33
right_hab_neuron_locs = findall(right_hab_neurons)
nRHabNeurons = minimum([nNeuronsForStim, size(right_hab_neuron_locs,1)])
right_hab_neuron_locs = sample(right_hab_neuron_locs, nRHabNeurons, replace=false)
left_hab_neuron_locs = findall(left_hab_neurons)
nLHabNeurons = minimum([nNeuronsForStim, size(left_hab_neuron_locs,1)])
left_hab_neuron_locs = sample(left_hab_neuron_locs, nLHabNeurons, replace=false)
raphe_neuron_locs = findall(raphe_neurons)
nRapheNeurons = minimum([nNeuronsForStim, size(raphe_neuron_locs,1)])
raphe_neuron_locs = sample(raphe_neuron_locs, nRapheNeurons, replace=false)

# TODO:
# - how to convert right_hab_neuron_locs to Array...?
# - transform into targets array, and call create_slm_stim
# - allow targets to pass 3D array with z-offset
# - create the file that Sean's randomized trials button makes (this way we control)
habZoffset = (zero_plane - planeHab) * zseries_step_size * 10^-6
right_hab_targets = vcat(Float64.(hcat(collect.(Tuple.(right_hab_neuron_locs))...)),
ones(Float64,1,nRHabNeurons)*habZoffset)'
# must be in 512x512 image space!! -> cropping makes life hard :/
right_hab_targets[:,[1,2]] ./= 2
right_hab_targets = right_hab_targets[:,[2,1,3]]
@warn "hack solution for the flip problem..."
right_hab_targets[:,1] .= 512 .- right_hab_targets[:,1]
right_hab_targets[:,2] .= 512 .- right_hab_targets[:,2]

left_hab_targets = vcat(Float64.(hcat(collect.(Tuple.(left_hab_neuron_locs))...)),
ones(Float64,1,nLHabNeurons)*habZoffset)'
left_hab_targets[:,[1,2]] ./= 2
left_hab_targets = left_hab_targets[:,[2,1,3]]
left_hab_targets[:,1] .= 512 .- left_hab_targets[:,1]
left_hab_targets[:,2] .= 512 .- left_hab_targets[:,2]

rapheZoffset = (zero_plane - planeRaphe) * 2. * 10^-6
raphe_targets = vcat(Float64.(hcat(collect.(Tuple.(raphe_neuron_locs))...)),
ones(Float64,1,nRapheNeurons)*rapheZoffset)'
raphe_targets[:,[1,2]] ./= 2
raphe_targets = raphe_targets[:,[2,1,3]]
raphe_targets[:,1] .= 512 .- raphe_targets[:,1]
raphe_targets[:,2] .= 512 .- raphe_targets[:,2]

create_slm_stim([left_hab_targets, right_hab_targets, raphe_targets],
    # "diff_num_neurons")
    # "fish2_with_resting_stim_lrhab_raphe")
    "fish2_with_resting_stim_lrhab_raphe_no_adj_h2b")
    # "fish2_with_resting_stim_lrhab_raphe_SyN")

## Visualize stim targets

im = RGB.(adj_zseries)
stim_points = zeros(Bool,size(adj_zseries))
stim_points[right_hab_neuron_locs,planeHab] .= true
stim_points[left_hab_neuron_locs,planeHab] .= true
stim_points[raphe_neuron_locs,planeRaphe] .= true
stim_points = dilate(dilate(stim_points))
channelview(im)[[1,3],:,:,:] .= 0
channelview(im)[1,:,:,:] .= float(stim_points)
channelview(im)[3,:,:,:] .= 0.5 * (hab_mask .| raphe_mask)
imshow(im)
@show planeHab, planeRaphe, zero_plane
println("Set first plane to  $(numZplanes-planeHab-1) and step size $step_size μm to image 3 planes")