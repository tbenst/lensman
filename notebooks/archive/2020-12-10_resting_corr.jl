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


##


# restingPreDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-pre-043"
# restingMidDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-mid-046/"
# restingPostDir = "/mnt/deissero/users/tyler/b115/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-resting-post-048"
restingPreDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_resting_pre-042"


tifDir = restingMidDir
expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
tseries = loadTseries(tifDir)
(H, W, Z, T) = size(tseries)


## noise correlation
cmax = 0.2
lag = 0
for stimNum in 1:3
    imageCorr = imageCorrWithMask(tseries[:,:,:,:],
        regionMasks[:,:,:,stimNum], lag=lag);
    # h5write(joinpath(fishDir,expName, "_imageCorr.h5"), "/imageCorr", imageCorr)
    fig = plt.figure(figsize=(50,10))
    plt.axis("off")
    plt.imshow(hcat([imageCorr[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
        clim=(-cmax, cmax))
    fig.savefig(joinpath(fishDir, "plots",expName*"_midCorr_$(lag)idx-lag_region$stimNum.png"))
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