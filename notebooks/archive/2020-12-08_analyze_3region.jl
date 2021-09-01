# need X11 forwarding for ssh, `ssh -X osprey`, then `echo $DISPLAY`, and fill below
# ENV["DISPLAY"] = "localhost:11"
ENV["DISPLAY"] = ":0"
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, Plots, HDF5,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman,
    Zarr
using Base.Iterators: peel
import Unitful: μm
##
# can set in Julia extension settings for vscode, search "num threads"
@assert Threads.nthreads() == 16

# can switch to this local NVME drive once rsync finishes... 
# dataDir = "/scratch/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab_raphe_40trial-045/"
# for now use network drive (tyler has to auth once per day...)
dataDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab_raphe_40trial-045/"


# TODO: copy slmExpDir (see below) to, say,  `/scratch/SLM_files/`, so no dependency on network
slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"

# painfully slow to load 90GB over 1Gb connection, so switch to /scratch once rsync finishes
tseries = loadTseries(dataDir)
(H, W, Z, T) = size(tseries);
##
voltageFile = glob("*VoltageRecording*.csv", dataDir)[1]
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)

## single stim example
imshow(imadjustintensity(tseries[:,:,:,stimStartIdx[1]-1:stimEndIdx[1]+1]))

## read slm stim files
dataFolders = splitpath(dataDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)

# if imaging many planes, may need to read another xml file since didn't get all planes in first file
@assert length(etlVals) == Z
@info "assume imaging from top-down"

#
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, _ = getTrialOrder(slmExpDir, expDate)


nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
@assert nTrials == size(stimStartIdx,1) # check TTL start-times match 

target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

# Array of matrices indicating x, y z (or was it y, x, z..?)
groupLocs = mapTargetGroupsToPlane(target_groups, etlVals)

# cells = makeCellsDF(targetLocs, stimStartIdx, stimEndIdx, trialOrder)

## avg stim effec
pre = Int(ceil(1*frameRate))
post = Int(ceil(1*frameRate))

# warmup JIT
@time trialAverage(tseries, stimStartIdx[1:3], stimEndIdx[1:3], trialOrder[1:3],
    pre=pre, post=post);
@time avgStim = trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder;
    pre=pre, post=post);

## STA
# h5write(joinpath(fishDir,expName, "_avgStim.h5"), "/block1", avgStim)
imshow(avgStim[:,:,:,1,:]) # trial type 1 (left hab)
imshow(avgStim[:,:,:,2,:]) # right hab
imshow(avgStim[:,:,:,3,:]) # raphe

## construct region masks
@info "hardcoded microscope units / targetSize"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
targetSizePx = (7μm * 14.4/25) / microscope_units[1]

regionMasks = constructGroupMasks(groupLocs, H, W, Z, targetSizePx=targetSizePx);
# imshow(regionMasks)

## clear stim tseries...
tseries = nothing
GC.gc()
sleep(3)
GC.gc() #look at memory and verify freed...
sleep(3)
GC.gc() #look at memory and verify freed...
sleep(3)
@everywhere GC.gc() #look at memory and verify freed...
##
restingPreDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-pre-043"
expName = splitpath(restingPreDir)[end]
fishDir = joinpath(splitpath(restingPreDir)[1:end-1]...)
tseriesPre = loadTseries(restingPreDir)
(H, W, Z, T) = size(tseriesPre)


## noise correlation
#warmup 
#TODO: set JULIA_NUM_THREADS=`nproc`
# lagged is better? of like 15?
imageCorr = imageCorrWithMask(tseriesPre[:,:,:,:], regionMasks[:,:,:,1]);
# h5write(joinpath(fishDir,expName, "_imageCorr.h5"), "/imageCorr", imageCorr)

imshow(imageCorr)
