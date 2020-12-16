restingPreDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-resting-pre-037"
tifDir = restingPreDir

expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
tseries = loadTseries(tifDir)
(H, W, Z, T) = size(tseries)


dataFolders = splitpath(tifDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
volRate = frameRate / Z

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

@info "hardcoded microscope units / targetSize"
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
targetSizePx = (7μm * 14.4/25) / microscope_units[1]

# imshow(regionMasks)
regionMasks = constructGroupMasks(groupLocs, H, W, Z, targetSizePx=targetSizePx);



## noise correlation
cmax = 0.1
lag = 0
for stimNum in 1:3
    imageCorr = imageCorrWithMask(tseries[:,:,:,:],
        regionMasks[:,:,:,stimNum], lag=lag);
    # h5write(joinpath(fishDir,expName, "_imageCorr.h5"), "/imageCorr", imageCorr)
    global fig = plt.figure(figsize=(50,10))
    plt.axis("off")
    plt.imshow(hcat([imageCorr[:,:,z] for z in 1:Z]...), cmap="RdBu_r")
    plt.imshow(hcat([imageCorr[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
        clim=(-cmax, cmax))
    fig.savefig(joinpath(fishDir, "plots",expName*"_midCorr_$(lag)idx-lag_region$(stimNum)_cmax$(cmax).png"))
end
fig
## colorbar
fig = plt.figure(figsize=(4,6), dpi=300)
plt.imshow(vcat([imageCorr[:,:,z] for z in 1:3]...), cmap="RdBu_r",
    clim=(-cmax, cmax))
cbar = plt.colorbar()
cbar.ax.set_ylabel("correlation")
fig.savefig(joinpath(fishDir, "plots",expName*"_midCorr_colorbar.png"),
    bbox_inches="tight")
fig
