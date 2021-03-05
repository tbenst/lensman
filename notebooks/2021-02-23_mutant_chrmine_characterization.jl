## compare 
# ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
import Gadfly
using Unitful: μm, m, s, mW

##
# offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
offset = float(uconvert(m, 0μm)) / m # when using SLM2 since 2020-02-?
zOffset = offset * 1e6
tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
# tseriesRootDir = "/data/dlab/b115"
# tseriesRootDir = "/mnt/deissero/users/tyler/b115"
# tseriesDir = "/data/dlab/b115/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-1024cell-32concurrent-4freq-054"
# tseriesDir = "/data/dlab/b115/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-256cell-8concurrent-4freq-055"
# can't fit in memory on lensman, so use deis
tseriesDir = joinpath(tseriesRootDir, "2021-02-16_6f_h33r_f0_6dpf/fish2/TSeries-256cell-8concurrent-4freq-051")
# tseriesDir = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043")
tseriesDir = joinpath(tseriesRootDir, "2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2/TSeries-128cell-4concurrent-3power-skip7-044")
tseriesDir = "$tseriesRootDir/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043"
tseriesDir = "$tseriesRootDir/2021-01-26_rsChRmine_6f_7dpf/fish2/TSeries-32concurrent-256trial-2rep-4power-045" # only one stim...?
tseriesDir = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-10zplane-077"
tseriesDir = "$tseriesRootDir/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060"


# possibly compare to...
# 2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/ (4power)
# 2021-02-15_wt_chrmine_gc6f/fish1/TSeries-1024cell-4freq-skip-first-066 (4freq; too large for memory on lensman)

# tyh5Path = tseriesDir

if occursin("freq", tseriesDir)
    exp_param = :stimFreq
elseif occursin("power", tseriesDir)
    exp_param = :laserPower
end

if tseriesDir[end] == "/"
    tseriesDir = tseriesDir[1:end-1]
end

fishDir = joinpath(splitpath(tseriesDir)[1:end-1]...)
expName = splitpath(tseriesDir)[end]
# tylerSLMDir = joinpath(fishDir, "slm")
tylerSLMDir = fishDir

tseries = loadTseries(tseriesDir);
# tyh5

# fishDir = joinpath(splitpath(tyh5Path)[1:end-1]...)
# expName = replace(splitpath(tyh5Path)[end], ".ty.h5" => "")
# tseriesDir = joinpath(fishDir, expName)

# tseries = h5read(tyh5Path, "/imaging/raw")
# @assert size(tseries,4)==1
# tseries = permutedims(tseries, (2,1,3,4,5))
# tseries = tseries[:,:,:,1,:];

##
(H, W, Z, T) = size(tseries)
@show (H, W, Z, T)
# slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
slmDir = "/oak/stanford/groups/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment/"
plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

# TODO: wrap more of cell in function to reduce inline code..?
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]

# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
Z = size(etlVals,1)
volRate = frameRate / Z
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

## read power
@assert length(glob("*.txt", tylerSLMDir)) == 1 # if not, need to be careful to choose
slmTxtFile = glob("*.txt", tylerSLMDir)[1]
stimGroupDF = CSV.File(open(read, slmTxtFile), header=["filepath", "powerFraction"]) |> DataFrame
stimGroupDF = stimGroupDF[trialOrder,:]


getMatStimFreq(mat) = sum((~).(sum.(mat["cfg"]["exp"]["targets"]) .≈ 0.0))*5
getSLMnum(mat) = size(mat["cfg"]["exp"]["targets"][1]) == (0,0) ? 1 : 2

mat = matread.(findMatGroups(slmExpDir)[1])
slmNum = getSLMnum(mat)
##
@warn "hardcoded laser power"
firstTargetGroup = matread.(findMatGroups(slmExpDir)[1])
powerPerCell = firstTargetGroup["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"]
if expDate < Date("2021-02-24")
    slm1Power = 570mW
    slm2Power = 380mW
elseif expDate <= Date("2021-03-01")
    slm1Power = 450mW
    slm2Power = 375mW
else
    slm1Power = 377mW
    slm2Power = 305mW
end
if slmNum == 1
    slmpowerPerCell = slm1Power * powerPerCell / 1000
elseif slmNum == 2
    slmpowerPerCell = slm2Power * powerPerCell / 1000
end

stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
allWithin(diff(stimStartIdx),0.05)
# adjust for expected number of stimuli....
# @assert length(stimStartIdx) == 32
@show length(stimStartIdx)
allWithin(stimEndIdx .- stimStartIdx,0.05)
avgStimDuration = mean(stimEndIdx .- stimStartIdx)
# adjust for volume rate / stim duration...
@assert 30 > avgStimDuration > 3

# Assumes no sequence of stim
# for 5Hz clock


target_groups = []
group_stim_freq = []
for mat in matread.(findMatGroups(slmExpDir))
    push!(target_groups, mat["cfg"]["maskS"]["targets"][1])
    push!(group_stim_freq, getMatStimFreq(mat))
end

# target_groups = [mat["cfg"]["maskS"]["targets"][1]
#     for mat in matread.(findMatGroups(slmExpDir))]

# group_stim_freq = [getMatStimFreq(mat) for mat in matread.(findMatGroups(slmExpDir))]

## TODO: very memory inefficient!  eats 90GB of RAM beyond tseries?!?
nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
# @assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}
is1024 = size(tseries,1)==1024
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
is1024=is1024, zOffset=zOffset)

# @warn "only using first 32... (also why only 63 stimuli??"
# cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx[1:32], stimEndIdx[1:32], trialOrder[1:32])
# cells2 = makeCellsDF(targetsWithPlaneIndex, stimStartIdx[33:63], stimEndIdx[33:63], trialOrder[33:63])
# trialOrder = trialOrder[7:end]
cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
cells[!, :stimFreq] = map(g->group_stim_freq[trialOrder][g], cells.stimGroup)
cells[!, :laserPower] = round.(typeof(1.0mW), map(g->stimGroupDF.powerFraction[g], cells.stimGroup) .* slmpowerPerCell, digits=1)

stimLocs = map(CartesianIndex ∘ Tuple, eachrow(cells[:,[2,1]]))

stimPlane = Int(mean(map(x->mean(x[:,3]), targetsWithPlaneIndex)))
avgImage = dropdims(mean(tseries[:,:,stimPlane,:], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
avgImageAdj = RGB.(avgImageAdj)
channelview(avgImageAdj)[[1,3],:,:,:] .= 0

if size(tseries,1)==512
    microscope_units = (2* 0.6299544139175637μm, 2* 0.6299544139175637μm, 2.0μm)
elseif size(tseries,1)==1024
    microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
else
    @error "unknown microscope units"
end

# assume that "7um" spiral galvo signal is calibrated for 16x, not 25x
# "16x" is actually 14.4x
targetSizePx = 7μm * (14.4/25) / microscope_units[1]
avgImgWithTargets = addTargetsToImage(copy(avgImageAdj), cartIdx2Array(stimLocs),
    targetSize=targetSizePx)
save(joinpath(plotDir,"$(expName)_avgImgWithTargets.png"), avgImgWithTargets)

winSize = Int(ceil(5*volRate))
# delay=Int(ceil(volRate*2))
delay=0

roiMask = constructROImasks(cells, H, W, Z, targetSizePx=targetSizePx);

cellsDF = makeCellsDF(tseries, cells, roiMask, winSize=winSize, delay=delay);
# cellsDF2 = makeCellsDF(tseries, cells2, roiMask, winSize=winSize, delay=delay);

rankings = sortperm(cellsDF.df_f, lt=(>))


meanDf_f = combine(groupby(cellsDF, :cellID), :df_f => mean).df_f_mean
cellIDrankings = sortperm(meanDf_f, lt=(>))

# rankings2 = sortperm(cellsDF2.df_f, lt=(>))
numNaN = sum(isnan.(cellsDF.df_f))
@assert numNaN == 0

@show sum(cellsDF.df_f.>0.1), sum(cellsDF.df_f.>0.4), sum(cellsDF.df_f.>1.)
# @show sum(cellsDF2.df_f.>0.1), sum(cellsDF2.df_f.>0.4), sum(cellsDF2.df_f.>1.)

# only 26 are >40% for two stim in control...
# intersect(cellsDF.x[cellsDF.df_f.>0.4], cellsDF2.x[cellsDF2.df_f.>0.4])
# @show first(cellsDF[rankings[numNaN+1:end],:],50)

open(tseriesDir*"_cellsDF.arrow", "w") do io
    tempDF = copy(cellsDF)
    tempDF[!, :laserPower] ./= 1mW
    Arrow.write(io, tempDF)
end



##
histogram(cellsDF.df_f)

##
# volRate = 30
before = Int(ceil(volRate*10))
after = Int(ceil(volRate*20))
# best
plots = []
nplots = 40

# nplots = 8
# for idx in rankings[end:-1:end-nplots+1]
# for idx in 1:nplots
@showprogress for cellID in cellIDrankings[1:nplots]
    indices = findall(cellsDF.cellID .== cellID)
    push!(plots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after, colorBy=exp_param))
end
p = plot(plots..., layout=(8,5), legend=true, size=(1024*4,1024*2))
savefig(p, joinpath(plotDir,"$(expName)_top40_traces_rollmedian.png"))
savefig(p, joinpath(plotDir,"$(expName)_top40_traces_rollmedian.svg"))
p
##

# TODO: why do we have NaN values...?
cmax = maximum(abs.(cellsDF.df_f[(~).(isnan.(cellsDF.df_f))]))
p_df_map = Gadfly.with_theme(:dark) do
    Gadfly.plot(cellsDF[(~).(isnan.(cellsDF.df_f)), :], x=:x, y=:y, Gadfly.Geom.point, size=[targetSizePx], color=:df_f,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true, fixed=true))
end

# img = PNG(joinpath(plotDir, "$(expName)_df_f_map.png"), 6inch, 5inch)
img = SVG(joinpath(plotDir, "$(expName)_df_f_map.svg"), 6inch, 5inch)
Gadfly.draw(img, p_df_map)
p_df_map

## extract and save traces for easy plotting later...
fluor = DataFrame(time=Float64[], f=Float64[], cellID=UInt32[], stimStart=UInt32[],
    laserPower=UInt32[], stimFreq=Float64[])

# TODO also compare on target vs off-target traces...? Number of off-target change
# as function of # of cell stim?

nTime = before + maximum(cells.stimStop - cells.stimStart) + after

for cell in eachrow(cells)
    x, y, z, stimStart, stimStop, laserPower, stimFreq, cellID = cell[[:x, :y, :z, :stimStart, :stimStop, :laserPower, :stimFreq, :cellID]]
    roiM = roiMask[cellID]
    theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
    theStart = maximum([stimStart-before, 1])
    plotRange = theStart:theEnd
    timeRange = (plotRange .- stimStart) ./ volRate
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    for (t,f) in zip(timeRange, fluorescentTrace)
        push!(fluor, (time=t, f=f, cellID=cellID, stimStart=stimStart,
            stimFreq=stimFreq, laserPower=laserPower))
    end
end
first(fluor,5)

open(tseriesDir*"_cells_fluorescence.arrow", "w") do io
    Arrow.write(io, fluor)
end


## 





## TODO: use 8vs32cell concurrent notebook as boilerplate to make WT vs rs vs H33R vs control plot...
# 8cell vs 32cell concurrent!