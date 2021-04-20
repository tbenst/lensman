# written for the 3 or 4 region 5-block stims
using FileIO, NRRD, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, Images, #ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, HDF5, Arrow, ProgressMeter,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman,
    NPZ
using ImageView
import Base.Threads.@threads
using Gadfly
import PyPlot
plt = PyPlot
np = pyimport("numpy")
sio = pyimport("scipy.io")

matplotlib = plt.matplotlib
##
zbrain_dir = "/data/zbrain"

# slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
slmDir = "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/"
# slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"
exp_dir = "2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1"
# fish_dir = "/data/dlab/b115/$exp_dir/"
fish_dir = "/mnt/deissero/users/tyler/b115/$exp_dir/"
tseries_dirs = glob("TSeries*", fish_dir)
# sort by recording number
tseries_dirs = tseries_dirs[sortperm(map(x->x[1], match.(r".*(\d{3})", tseries_dirs)))]
# when suite2p processes multiple files, output in first folder
first_expName = splitpath(tseries_dirs[1])[end]
suite2p_dir = "/data/dlab/b115/process-output/$exp_dir/$first_expName/hdf5/data/suite2p/"
hdf5_dir = "/data/dlab/b115/process-output/$exp_dir/$first_expName/hdf5/"
# read timesteps per experiment block
blocks_timesteps = zeros(Int64, length(tseries_dirs))
global Z = 0
for (i,tseries_dir) in enumerate(tseries_dirs)
    _,_,Z,T,_ = tseriesTiffDirMetadata(tseries_dir)
    blocks_timesteps[i] = T
end

## read nCells
plane_dirs = glob("plane*", suite2p_dir)

iscell = npzread(joinpath(suite2p_dir, "combined", "iscell.npy"));
nCells = size(iscell,1)

combined_stat = PyObject[]
for z in 1:Z
    # have to process suite2p plane by plane due to
    # https://github.com/MouseLand/suite2p/issues/655
    stat = np.load(joinpath(suite2p_dir, "plane$(z-1)", "stat.npy"), allow_pickle=true);
    for k in keys(stat)
        # uses 1 indexing (julia) rather than 0 indexing (suite2p/python)
        stat[k].update(iplane=z) # can't do []  assignment in PyCall
    end
    cell_idx_offset = length(combined_stat)
    combined_stat = vcat(combined_stat, stat)
end
cell_centers = map(x->np.array(get(x,"med")), combined_stat);


##
@info "assumes 2nd & 4th are stim blocks"
# get cell masks per z-plane

# output:
# - stim: nCell x nStim array, with power (mW) per cell
# - stim_start: nStim array, with timestep as value
combined_stim_start = Int64[]
combined_stim_end = Int64[]
combined_trial_order = Int64[]
cell_in_mask = zeros(Bool,nCells,nStimuli)
prev_stim_masks = nothing
for (i,tseries_dir) in enumerate(tseries_dirs[[2,4]])
    # get stim start times & add offset
    if i == 1
        t_offset = blocks_timesteps[1]
    elseif i == 2
        t_offset = sum(blocks_timesteps[1:3])
    else
        error("not implemented")
    end

    H,W,Z,T,_ = tseriesTiffDirMetadata(tseries_dir)
    expName = splitpath(tseries_dir)[end]
    voltageFile = joinpath(tseries_dir, expName*"_Cycle00001_VoltageRecording_001.csv")
    if ~isfile(voltageFile)
        # avoid scanning list of many files...
        voltageFile = glob("*VoltageRecording*.csv", tseries_dir)[1]
    end

    stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
    stimStartIdx .+= t_offset
    stimEndIdx .+= t_offset
    combined_stim_start = vcat(combined_stim_start, stimStartIdx)
    combined_stim_end = vcat(combined_stim_end, stimEndIdx)
    dataFolders = splitpath(tseries_dir)
    xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
    expDate, frameRate, etlVals = getExpData(xmlPath)
    volRate = frameRate / Z
    slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
    trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)
    combined_trial_order = vcat(combined_trial_order, trialOrder)
    nStimuli = maximum(trialOrder)
    nTrials = size(trialOrder,1)
    @assert nTrials == size(stimStartIdx,1) # check TTL start-times match 
    target_groups = [mat["cfg"]["maskS"]["targets"][1]
        for mat in matread.(findMatGroups(slmExpDir))]
    target_groups_with_plane = mapTargetGroupsToPlane(target_groups, etlVals,
        is1024=(W == 1024))
    target_groups_with_plane = map(x->Int.(round.(x)), target_groups_with_plane)
    stim_masks = constructGroupMasks(target_groups_with_plane, H, W, Z)
    if ~isnothing(prev_stim_masks)
        @assert all(stim_masks .== prev_stim_masks) "stim masks differ between blocks"
    end

    for (c,center) in enumerate(cell_centers)
        # uses 1 indexing (julia) rather than 0 indexing (suite2p/python)
        z = combined_stat[c].get("iplane")
        for s in nStimuli
            cell_in_mask[c,s] = stim_masks[center[[2,1]]...,z,s]
            # cell_in_mask[c,s] = stim_masks[center...,z,s]
        end     
    end
    prev_stim_masks = stim_masks
end
@show sum(cell_in_mask,dims=1)[1,:]

## TODO: save stim_masks and trial order and start times
stim_h5 = h5open(joinpath(hdf5_dir, "stim.h5"), "w")
stim_h5["stim_start_idx"] = combined_stim_start
stim_h5["stim_end_idx"] = combined_stim_end
stim_h5["trial_order"] = combined_trial_order
stim_h5["stim_masks"] = prev_stim_masks

## get zbrain masks
affine_transform_path = glob("*GenericAffine.mat", fish_dir)[1]
affine_transform = sio.loadmat(affine_transform_path)
warpedname = glob("*_Warped.nii.gz", fish_dir)[1]
zbrain_registered = niread(warpedname);

# adj_zbrain_registered = adjust_histogram(imadjustintensity(zbrain_registered), Equalization());
adj_zbrain_registered = adjust_gamma(imadjustintensity(zbrain_registered), 0.2);

## visualize results as sanity check
im = RGB.(adj_zbrain_registered)
channelview(im)[2:3,:,:,:] .= 0
# channelview(im)[2,:,:,:] .= adj_zseries
imshow(im)
masks = matread("$zbrain_dir/MaskDatabase.mat")

## TODO: for each mask, convert with antsApplyTransforms, and save a new HDF5 file for Zseries
# also need to record which planes correspond to our imaging planes

## debug viz
imshow(stim_masks)
cell_planes = map(x->get(x,"iplane"), combined_stat);
cell_groups = vcat([transpose([c[[2,1]]..., z]) for (c,z) in zip(cell_centers,cell_planes)]...);
cell_masks = constructGroupMasks([cell_groups], H, W, Z);
imshow(cell_masks)