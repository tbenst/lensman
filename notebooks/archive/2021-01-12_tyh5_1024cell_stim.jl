##
ENV["DISPLAY"] = "localhost:14.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter
import Gadfly
using Unitful: μm, m, s

##
# offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
offset = float(uconvert(m, 0μm)) / m # when using SLM2 since 2020-02-?
zOffset = offset * 1e6

# tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-047.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-047.ty.h5"
# tseriesDir = "/mnt/deissero/users/tyler/b115/2021-01-11_chrmine-kv2.1_h2b6s_6dpf/fish1_nochrmine/TSeries-1024cell-32concurrent-043"
# tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-18_chrmine_kv2.1_h2b6s_6dpf/fish1_chrmine/TSeries-1024cell-32concurrent-2reps-043.ty.h5"
# tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-18_chrmine_kv2.1_h2b6s_6dpf/fish2_nochrmine/TSeries-1024cell-32concurrent-2reps-044.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-5ppc-048.ty.h5"
# tyh5Path = "/mnt/deissero/users/tyler/b115/2020-12-17_gc6f-h33r_7dpf/fish1/TSeries_10stim_4times-40power-newFOV-050.ty.h5"
# tseriesDir = "/mnt/deissero/users/tyler/b115/2020-12-15_h2b6s_chrmine_kv2.1_5dpf/fish2/TSeries-1024cell-32concurrent-038"
# tseriesDir = "/mnt/deissero/users/tyler/b115/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-1024cell-32concurrent-4freq-054"
tseriesDir = "/data/dlab/b115/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-1024cell-32concurrent-4freq-054"
# tyh5Path = tseriesDir

if tseriesDir[end] == "/"
    tseriesDir = tseriesDir[1:end-1]
end


tseries = loadTseries(tseriesDir);
fishDir = joinpath(splitpath(tseriesDir)[1:end-1]...)
expName = splitpath(tseriesDir)[end]

# tyh5

# fishDir = joinpath(splitpath(tyh5Path)[1:end-1]...)
# expName = replace(splitpath(tyh5Path)[end], ".ty.h5" => "")
# tseriesDir = joinpath(fishDir, expName)

# tseries = h5read(tyh5Path, "/imaging/raw")
# @assert size(tseries,4)==1
# tseries = permutedims(tseries, (2,1,3,4,5))
# tseries = tseries[:,:,:,1,:];

## plot whole image...
brainTrace = mean(tseries, dims=(1,2,3))[1,1,1,:]
brainTrace = imageJkalmanFilter(brainTrace)
p = plot((1:size(brainTrace,1)) ./ volRate, brainTrace, left_margin=50px,
    xlabel="time (s)", ylabel="fluorescence")
for (stimStart,stimStop) in zip(stimStartIdx, stimEndIdx)
    plot!(Lensman.rectangle((stimStop-stimStart)/volRate,maximum(brainTrace)-minimum(brainTrace),
        stimStart/volRate, minimum(brainTrace)), opacity=.5, label="")
end
p

##
(H, W, Z, T) = size(tseries)
@show (H, W, Z, T)
slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
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

@warn "hardcoded laser power"
firstTargetGroup = matread.(findMatGroups(slmExpDir)[1])
powerPerCell = firstTargetGroup["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"]
# slm1Power = 850
# slm1powerPerCell = slm1Power * powerPerCell / 1000


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
getMatStimFreq(mat) = sum((~).(sum.(mat["cfg"]["exp"]["targets"]) .≈ 0.0))*5
target_groups = []
group_stim_freq = []
for mat in matread.(findMatGroups(slmExpDir))
    push!(target_groups, mat["cfg"]["maskS"]["targets"][1])
    push!(group_stim_freq, getMatStimFreq(mat))
end

# target_groups = [mat["cfg"]["maskS"]["targets"][1]
#     for mat in matread.(findMatGroups(slmExpDir))]

# group_stim_freq = [getMatStimFreq(mat) for mat in matread.(findMatGroups(slmExpDir))]

##
exMat = matread(findMatGroups(slmExpDir)[1])["cfg"]
println(keys(exMat))
exMat["exp"]

##
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
cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
cells[!, :stimFreq] = map(g->group_stim_freq[trialOrder][g], cells.stimGroup)

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
    Arrow.write(io, cellsDF)
end



##
histogram(cellsDF.df_f)

##
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
    push!(plots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after, colorBy=:stimFreq))
end
p = plot(plots..., layout=(8,5), legend=true, size=(1024*4,1024*2))
savefig(p, joinpath(plotDir,"$(expName)_top40_traces.png"))
savefig(p, joinpath(plotDir,"$(expName)_top40_traces.svg"))
p
##

# TODO: why do we have NaN values...?
cmax = maximum(abs.(cellsDF.df_f[(~).(isnan.(cellsDF.df_f))]))
p_df_map = Gadfly.with_theme(:dark) do
    Gadfly.plot(cellsDF[(~).(isnan.(cellsDF.df_f)), :], x=:x, y=:y, Gadfly.Geom.point, size=[targetSizePx], color=:df_f,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true, fixed=true))
end

# img = PNG(joinpath(plotDir, "$(expName)_df_f_map_below10.png"), 6inch, 5inch)
# Gadfly.draw(img, p_df_map)
p_df_map

## extract and save traces for easy plotting later...
fluor = DataFrame(time=Float64[], f=Float64[], cellID=UInt32[], stimStart=UInt32[])

# TODO also compare on target vs off-target traces...? Number of off-target change
# as function of # of cell stim?

nTime = before + maximum(cells.stimStop - cells.stimStart) + after

for cell in eachrow(cells)
    x, y, z, stimStart, stimStop, stimFreq, cellID = cell[[:x, :y, :z, :stimStart, :stimStop, :stimFreq, :cellID]]
    roiM = roiMask[cellID]
    theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
    theStart = maximum([stimStart-before, 1])
    plotRange = theStart:theEnd
    timeRange = (plotRange .- stimStart) ./ volRate
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    for (t,f) in zip(timeRange, fluorescentTrace)
        push!(fluor, (time=t, f=f, cellID=cellID, stimStart=stimStart))
    end
end
first(fluor,5)

open(tseriesDir*"_cells_fluorescence.arrow", "w") do io
    Arrow.write(io, fluor)
end