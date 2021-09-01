# MUST RUN ON LINUX FOR ANTS support (unless compiled for windows...)

# This notebook takes a 820nm Zseries, and:
# - registers to Zbrain
# - finds best plane for imaging raphe + hab
# - finds centers of neurons
# - saves this data in case session crashes
# - creates .mat file for use with Sean's multiSLM software stims

## run cell to init VSCode Julia
@warn raw"ensure using Olympus 25x in Prairie View, else ETL will be off"
##
# ENV["DISPLAY"] = "localhost:12.0"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates, DataFrames, CSV
import Base.Threads.@threads
using Unitful: μm, m, s
sio = pyimport("scipy.io")
if Sys.iswindows()
    error("this script only supports linux!")
end
zbrain_dir = "/data/zbrain"
offset = float(uconvert(m, 0μm)) / m # since 2020-02-10 for SLM2
slmNum = 2
zOffset = offset * 1e6

##
# fishDir = "/mnt/b115_data/tyler/2021-02-02_wt_chrmine_GC6f/fish4"
# fishDir = "/mnt/b115_data/tyler/2021-02-09_gcamp6f_7dpf/fish2"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-16_6f_h33r_f0_6dpf/fish2"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-04-13_wt-chrmine_6dpf_h2b6s/fish1"
# fishDir = "/scratch/b115/2021-05-04_hsChRmine_h2b6s_6dpf/fish1"
fishDir = "/scratch/b115/2021-05-05_wt-chrmine_6f_6dpf/fish1"

zseriesDir = glob("ZSeries*", fishDir)[end]

zseries = readZseriesTiffDir(zseriesDir)
H, W, Z = size(zseries)

lateral_unit = microscope_lateral_unit(W)
microscope_units = (lateral_unit, lateral_unit, 2.0μm)

H, W, Z
##
# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), microscope_units)
# adj_zseries = adjust_histogram(imadjustintensity(zseries), Equalization())
adj_zseries = imadjustintensity(adjust_gamma(zseries,0.5));
adj_zseries = AxisArray(adj_zseries, (:y, :x, :z), microscope_units);
@assert microscope_units[1] > 0.1μm
microscope_units

#####
## Read Zbrain cytosolic gcamp
#####
@warn "update for H2B vs cytosolic!!"
gcamp_zbrain = AxisArray(permutedims(
        # h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-GCaMP5G_6dpf_MeanImageOf7Fish"),
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), zbrain_units)
gcamp_zbrain = AxisArray(permutedims(
        load("$zbrain_dir/Ref20131120pt14pl2.nrrd"),
        (1,2,3)),
    (:x, :y, :z), zbrain_units);

gcamp_zbrain = reverse(gcamp_zbrain,dims=2); # face right
gcamp_zbrain = reverse(gcamp_zbrain,dims=3); # top to bottom

ps = pixelspacing(gcamp_zbrain); SD = ps./maximum(ps)
# adj_gcamp_zbrain = adjust_histogram(imadjustintensity(gcamp_zbrain), Equalization())
# adj_gcamp_zbrain = adjust_gamma(imadjustintensity(gcamp_zbrain), 0.2)
##

adj_gcamp_zbrain = adjust_gamma(gcamp_zbrain, 1.1);
# adj_gcamp_zbrain = adjust_gamma(gcamp_zbrain, 0.2);
adj_gcamp_zbrain = AxisArray(adj_gcamp_zbrain, (:y, :x, :z), zbrain_units);
size(gcamp_zbrain)
imshow(adj_gcamp_zbrain)
# ## Examine our registration starting point...
# im = RGB.(adjust_histogram(imadjustintensity(gcamp_zbrain), Equalization()))
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
# fixed = zseries
fixed = adj_zseries
moving = gcamp_zbrain
# register_method = "SyN" # 1000 seconds, no improvement
register_method = "affine" # 150 seconds
# moving = adj_gcamp_zbrain
fixedname = joinpath(fishDir, "zseries_ours.nrrd")
save(fixedname, fixed)
movingname = joinpath(fishDir, "zbrain_atlas.nrrd")
save(movingname, moving)
outprefix = Dates.format(DateTime(now()), DateFormat("YmmddTHHMMSSsZ"))
if register_method == "SyN"
    outname = joinpath(fishDir, outprefix*"_SyN")
else
    outname = joinpath(fishDir, outprefix)
end
println("creating files with name: $outname")
# no SyN so faster, but kinda meh registration...

# originally was -r instead of -q
# maybe try both?
# -q is --initial-fixed-transform
# -r is --initial-moving-transform

# AFFINE ONLY
# from ANTs zfish atlas paper
# interpolation = "WelchWindowedSinc"
interpolation = "BSpline"
histmatch = 0
# winsorize = ""
# initial_moving_type = 0
# sampling_frac = 0.25
maxiter = 200
threshold = "1e-6"
initial_moving_type = 1
sampling_frac = 1.0
synThreshold = "1e-6"
if register_method == "SyN"
    syn_cmd = `-t SyN\[0.05,6,0.5\] -m CC\[$fixedname, $movingname,1,2\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$synThreshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox`
else
    syn_cmd = ``
end
cmd = `/opt/ANTs/install/bin/antsRegistration -d 3 --float -o \[$outname\_, $outname\_Warped.nii.gz, $outname\_WarpedInv.nii.gz\] --interpolation $interpolation --use-histogram-matching $histmatch -r \[$fixedname, $movingname,$initial_moving_type\] -t rigid\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,    $sampling_frac\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox -t Affine\[0.1\] -m MI\[$fixedname, $movingname,1,32, Regular,$sampling_frac\] -c \[$maxiter\x$maxiter\x$maxiter\x0,$threshold,10\] --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox $syn_cmd`
println("calling ANTs...")
@time println(read(cmd, String))
##
# affine_transform_path = glob(outprefix*"_SyN"*"*GenericAffine.mat", fishDir)[1]
affine_transform_path = glob(outprefix*"*0GenericAffine.mat", fishDir)[1]
affine_transform = sio.loadmat(affine_transform_path)
if register_method=="SyN"
    syn_transform_path = glob(outprefix*"*1Warp.nii.gz", fishDir)[1]
    syn_transform = niread(syn_transform_path)
end
# warpedname = joinpath(outname*".nrrd")
warpedname = joinpath(outname*"_Warped.nii.gz")

# if this fails, tmpfs may be out of space
# try `rm -rf /run/user/1000/tyler/ANTs/`
zbrain_registered = niread(warpedname);
adj_zbrain_registered = adjust_histogram(imadjustintensity(zbrain_registered), Equalization());
println("finished registering!")

## visualize results as sanity check
im = RGB.(adj_zbrain_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
imshow(im)

## Read mask & transform to Zseries
# TODO: define our own hemisphere masks (split zbrain down the middle..?)
masks = matread("$zbrain_dir/MaskDatabase.mat")
raphe_mask = getMaskNameIncludes(masks, "Raphe");
raphe_mask = reverse(raphe_mask,dims=2); # face right
raphe_mask = reverse(raphe_mask,dims=3); # top to bottom
# make raphe binary again (weird annoying error but can run in shell...?)
raphe_mask = antsApplyTransforms2(fixedname, raphe_mask, affine_transform_path) .> 0;

hab_mask = getMaskNameIncludes(masks, "Habenula");
hab_mask = reverse(hab_mask,dims=2); # face right
hab_mask = reverse(hab_mask,dims=3); # top to bottom

hab_mask = antsApplyTransforms2(fixedname, hab_mask, affine_transform_path) .> 0;

## neuron locations by plane, takes ~30 seconds
neuronCenterMask = zeros(Bool, size(zseries)...)
for z in 1:size(neuronCenterMask, 3)
    neuronCenterMask[:,:,z] = findNeurons(zseries[:,:,z])
end

raphe_nNeurons_by_plane = mapslices(sum, raphe_mask .& neuronCenterMask, dims=[1,2])
max_raphe_nNeurons_in_plane = maximum(raphe_nNeurons_by_plane)
hab_nNeurons_by_plane = mapslices(sum, hab_mask .& neuronCenterMask, dims=[1,2])
max_hab_nNeurons_in_plane = maximum(hab_nNeurons_by_plane)
# auto determine plane
planeHab = Tuple(argmax(hab_nNeurons_by_plane))[3]
planeRaphe = Tuple(argmax(raphe_nNeurons_by_plane))[3]

# @warn "manual plane selection"
# planeHab = 104
# planeRaphe = 
##

fracPlaneZero = planeRaphe - (planeRaphe - planeHab)*.7
planeZero = Int(round(fracPlaneZero))
# -40 to 60 is acceptable range...
# +1 to move down one for spectral offset
numZplanes = size(zseries,3)
zseries_step_size = 2
zEtlZero = (fracPlaneZero + 1) * zseries_step_size
zHab = (planeHab + 1)*zseries_step_size
zRaphe = (planeRaphe + 1) * zseries_step_size
etlHab = zHab - zEtlZero # negative is up 
seanHab = etlHab + zOffset # (for sean's code)
etlRaphe = zRaphe - zEtlZero # positive is down 
seanRaphe = etlRaphe + zOffset

@info "assuming each plane is 2μm apart"
numImagingPlanes = 10
step = (etlRaphe - etlHab)/(numImagingPlanes-1)

step_size = step * zseries_step_size

# note that our zseries is inverted from PV...
# TODO: flip this so more intuitive...?
# TODO: why are we flipped lr...?
# println("Set Z 'zero' to be plane $(numZplanes-planeZero) and step size $step_size μm to image 3 planes")
@info "Ensure that set to to Olympus-25x, or else ETL will be wrong!"
# TODO: read state from xml file and communicate with PV so we get zseries microns correct...?

# MAKE PROTOCOL
# we image 3 planes: best hab, best raphe, and best compromise
# we have three stim groups:
#   (1) left hab (best hab plane) 
#   (2) right hab (best hab plane) 
#   (3) raphe (best raphe plane) 
plane_neurons = neuronCenterMask[:,:,planeHab]
right_mask = zeros(Bool, size(plane_neurons)...)
dividing_line = Int(round(mean(mask2ind(raphe_mask)[:,1])))
@info "dividing_line set to $dividing_line, check if reasonable?"
right_mask[1:dividing_line,:] .= 1
left_mask = zeros(Bool, size(plane_neurons)...)
left_mask[dividing_line+1:end,:] .= 1
hab_neurons = (plane_neurons .& hab_mask[:,:,planeHab])
right_hab_neurons = hab_neurons .& right_mask
left_hab_neurons = hab_neurons .& left_mask
@show sum(right_hab_neurons)
@show sum(left_hab_neurons)

plane_neurons = neuronCenterMask[:,:,planeRaphe]
raphe_neurons = plane_neurons .& raphe_mask[:,:,planeRaphe]
@show sum(raphe_neurons)

# sample consistent number of neurons
nNeuronsForStim = 50
right_hab_neuron_locs = findall(right_hab_neurons)
nRHabNeurons = minimum([nNeuronsForStim, size(right_hab_neuron_locs,1)])
right_hab_neuron_locs = sample(right_hab_neuron_locs, nRHabNeurons, replace=false)

left_hab_neuron_locs = findall(left_hab_neurons)
nLHabNeurons = minimum([nNeuronsForStim, size(left_hab_neuron_locs,1)])
left_hab_neuron_locs = sample(left_hab_neuron_locs, nLHabNeurons, replace=false)

raphe_neuron_locs = findall(raphe_neurons)
nRapheNeurons = minimum([nNeuronsForStim, size(raphe_neuron_locs,1)])
raphe_neuron_locs = sample(raphe_neuron_locs, nRapheNeurons, replace=false)

@warn "hardcode control region"
botLeftX = 200
botLeftY = 925
control_neuron_locs = CartesianIndex{2}[]
for i in 1:nNeuronsForStim
    # push!(control_neuron_locs, CartesianIndex(rand(botLeftY:1024), rand(1:botLeftX)))
    push!(control_neuron_locs, CartesianIndex(rand(50:150), rand(50:150)))
end
control_neuron_locs

# TODO:
# - how to convert right_hab_neuron_locs to Array...?
# - transform into targets array, and call create_slm_stim
# - allow targets to pass 3D array with z-offset
# - create the file that Sean's randomized trials button makes (this way we control)
habZoffset = seanHab * 10^-6
right_hab_targets = vcat(Float64.(hcat(collect.(Tuple.(right_hab_neuron_locs))...)),
ones(Float64,1,nRHabNeurons)*habZoffset)'
# must be in 512x512 image space!! -> cropping makes life hard :/
right_hab_targets[:,[1,2]] ./= 2
right_hab_targets = right_hab_targets[:,[2,1,3]]
# @info "hack solution for the flip problem..."
# right_hab_targets[:,1] .= 512 .- right_hab_targets[:,1]
# right_hab_targets[:,2] .= 512 .- right_hab_targets[:,2]

left_hab_targets = vcat(Float64.(hcat(collect.(Tuple.(left_hab_neuron_locs))...)),
ones(Float64,1,nLHabNeurons)*habZoffset)'
left_hab_targets[:,[1,2]] ./= 2
left_hab_targets = left_hab_targets[:,[2,1,3]]
# left_hab_targets[:,1] .= 512 .- left_hab_targets[:,1]
# left_hab_targets[:,2] .= 512 .- left_hab_targets[:,2]

rapheZoffset = seanRaphe * 10^-6
raphe_targets = vcat(Float64.(hcat(collect.(Tuple.(raphe_neuron_locs))...)),
ones(Float64,1,nRapheNeurons)*rapheZoffset)'
raphe_targets[:,[1,2]] ./= 2
raphe_targets = raphe_targets[:,[2,1,3]]
# raphe_targets[:,1] .= 512 .- raphe_targets[:,1]
# raphe_targets[:,2] .= 512 .- raphe_targets[:,2];

control_targets = vcat(Float64.(hcat(collect.(Tuple.(control_neuron_locs))...)),
ones(Float64,1,nLHabNeurons)*habZoffset)'
control_targets[:,[1,2]] ./= 2
control_targets = control_targets[:,[2,1,3]]


# Visualize stim targets
im = RGB.(adj_zseries)
stim_points = zeros(Bool,size(adj_zseries))
stim_points[right_hab_neuron_locs,planeHab] .= true
stim_points[left_hab_neuron_locs,planeHab] .= true
stim_points[raphe_neuron_locs,planeRaphe] .= true
stim_points = dilate(dilate(stim_points))
channelview(im)[[1,3],:,:,:] .= 0
channelview(im)[1,:,:,:] .= float(stim_points)
channelview(im)[3,:,:,:] .= 0.5 * (hab_mask .| raphe_mask)
# imshow(im)
@show planeHab, planeRaphe
# we move down one plane (2um) to account for spectral shift
println("Move Z to $(zEtlZero)μm")
println("Move ETL to hab plane @ $(etlHab)μm) & set as first")
println("Set step size to $(step)μm")
println("Set end at Raphe plane @ $(etlRaphe)μm")
#
step5 = Int(round(Int,step/2))
implanes = [planeHab, planeZero, planeRaphe]
imshow(im[:,:,implanes])

# zeroplane: Z=118
# ETL needs to go down 46 to reach raphe plane, not step5*2 = 32
# ETL needs to go down ~12 to reach zero plane
# ETL needs to go up -26 to reach hab plane, not -step5*2 = -32

# zeroplane should actually be 118-14 = 

# planeHab: Z=90
# planeZero: Z=130
# planeRaphe: Z=160

# 46: 90
# 62: 122
# 78: 154
##
magicX, magicY = markpoints_magic_numbers("B115")
write_markpoints(left_hab_neuron_locs, "$(outname)_left_hab.gpl",
    W=W, magicX=magicX, magicY=magicY, spiral_size=0.09)
write_markpoints(right_hab_neuron_locs, "$(outname)_right_hab.gpl",
    W=W, magicX=magicX, magicY=magicY, spiral_size=0.09)
write_markpoints(raphe_neuron_locs, "$(outname)_raphe.gpl",
    W=W, magicX=magicX, magicY=magicY, spiral_size=0.09)
write_markpoints(control_neuron_locs, "$(outname)_control.gpl",
    W=W, magicX=magicX, magicY=magicY, spiral_size=0.09)

##
picked_left_hab_neuron_locs = read_gpl(joinpath(fishDir, "handpicked_left_hab.gpl"),
    width=1024, height=1024)
picked_right_hab_neuron_locs = read_gpl(joinpath(fishDir, "handpicked_right_hab.gpl"),
    width=1024, height=1024)
picked_raphe_neuron_locs = read_gpl(joinpath(fishDir, "handpicked_raphe.gpl"),
    width=1024, height=1024)

right_hab_targets = vcat(Float64.(hcat(collect.(Tuple.(picked_right_hab_neuron_locs))...)),
ones(Float64,1,nRHabNeurons)*habZoffset)'
# must be in 512x512 image space!! -> cropping makes life hard :/
right_hab_targets[:,[1,2]] ./= 2
right_hab_targets = right_hab_targets[:,[2,1,3]]
# @info "hack solution for the flip problem..."
# right_hab_targets[:,1] .= 512 .- right_hab_targets[:,1]
# right_hab_targets[:,2] .= 512 .- right_hab_targets[:,2]

left_hab_targets = vcat(Float64.(hcat(collect.(Tuple.(picked_left_hab_neuron_locs))...)),
ones(Float64,1,nLHabNeurons)*habZoffset)'
left_hab_targets[:,[1,2]] ./= 2
left_hab_targets = left_hab_targets[:,[2,1,3]]
# left_hab_targets[:,1] .= 512 .- left_hab_targets[:,1]
# left_hab_targets[:,2] .= 512 .- left_hab_targets[:,2]

rapheZoffset = seanRaphe * 10^-6
raphe_targets = vcat(Float64.(hcat(collect.(Tuple.(picked_raphe_neuron_locs))...)),
ones(Float64,1,nRapheNeurons)*rapheZoffset)'
raphe_targets[:,[1,2]] ./= 2
raphe_targets = raphe_targets[:,[2,1,3]]


# Visualize stim targets
im = RGB.(adj_zseries)
stim_points = zeros(Bool,size(adj_zseries))
stim_points[picked_right_hab_neuron_locs,planeHab] .= true
stim_points[picked_left_hab_neuron_locs,planeHab] .= true
stim_points[control_neuron_locs,planeHab] .= true
stim_points[picked_raphe_neuron_locs,planeRaphe] .= true
stim_points = dilate(dilate(stim_points))
channelview(im)[[1,3],:,:,:] .= 0
channelview(im)[1,:,:,:] .= float(stim_points)
channelview(im)[3,:,:,:] .= 0.5 * (hab_mask .| raphe_mask)


imshow(im[:,:,implanes])
## SAVE

slmOutName = "$(outname)_lrhab_raphe_control"
if isfile(slmOutName*".txt")
    error("file already exists! refusing to clobber")
end

if fishDir[1:13] == "/mnt/deissero"
    localToRemote = matpath -> "O:\\" * replace(matpath[14:end], "/" => "\\")
elseif fishDir[1:12] == "/mnt/b115_data"
    localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\")
elseif fishDir[1:13] == "/scratch/b115"
    localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\")
else
    error("need to manually specify remote path.")
end
    
create_slm_stim([left_hab_targets, right_hab_targets, raphe_targets,
    control_targets],    slmOutName,
    localToRemote = localToRemote,
    powers=[1], slmNum=slmNum)

## trialOrder
nStims = 4
nReps = 8 # transition reps
nTrials = nStims^2 * nReps + 1

transitionsLeft = ones(4,4) * nReps
trialOrder = Int64[]

push!(trialOrder, rand(1:nStims))
for i in 2:nTrials
    prevStim = trialOrder[i-1]
    candidateTransitionsLeft = findall(transitionsLeft[prevStim,:] .> 0)
    nextStim = rand(candidateTransitionsLeft)
    push!(trialOrder, nextStim)
    transitionsLeft[nextStim,prevStim] -= 1
end
@assert length(trialOrder) == nTrials
##
trialOrderDF = DataFrame(copy(hcat(collect(1:nTrials), trialOrder)'))

CSV.write(outname*"_trialOrder.txt", trialOrderDF, header=false, delim="\t")
println("wrote $(outname*"_trialOrder.txt")")
println("be sure to modify mSLM/SetupFiles/Experiments/<TODAY>/trialOrder.txt")
