ENV["DISPLAY"] = "localhost:11"

using ImageView, Lensman, Glob, Dates, MAT
import Plots
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
##
tifDir = "/mnt/deissero/users/tyler/b115/2020-12-17_h2b6s_chrmine-kv2.1_7dpf/fish2/TSeries_cstoner_n64_b2_r8-046"
slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
tseries = loadTseries(tifDir);
(H, W, Z, T) = size(tseries)
##
voltageFile = glob("*VoltageRecording*.csv", tifDir)[1]
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
@show length(stimStartIdx)

dataFolders = splitpath(tifDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
volRate = frameRate / Z
@assert length(etlVals) == Z

slmExpDateDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDateDir, expDate)
@show nStimuli = maximum(trialOrder)
@show nTrials = size(trialOrder,1)
@show nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)

@assert length(findMatGroups(slmExpDir)) == nStimuli

target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals)

cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
xyzToIdx = Dict( (x,y,z) => i for (i,(x,y,z)) in
    enumerate(eachrow(unique(cells[!,[:x,:y,:z]]))) )
cells = addCellIdx(cells, xyzToIdx)
targetCellIdxs = map(targets->
    map(xyz->xyzToIdx[(xyz...)], eachrow(targets)),
    targetsWithPlaneIndex)

concurrencyMat = count_concurrency(targetCellIdxs)
# need to put eyes on to see if it worked...?
entangledImaps = entangledInfluenceMaps(cells, tseries);
targetCellIdxs
##
# 1. calculate df/f for each stimulus--the influence map of the group
# 2. construct concurrency matrix (X)
# 4. create entangled influence map for each stimulated neuron (Y)
# 5. minimize sum.(abs.(A)) subject to Y = AX 

pre = Int(ceil(2*volRate))
post = Int(ceil(2*volRate))

##  

transform(cells, [:x,:y,:z] => (xyz -> xyzToIdx[xyz]))
transform(cells, :x => (x -> x/512))
@show length(unique(eachrow(vcat(targetsWithPlaneIndex...))))