##
ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall
import Gadfly
using Unitful: μm, m, s, mW

##
offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
zOffset = offset * 1e6

# tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-18_chrmine_kv2.1_h2b6s_6dpf/fish1_chrmine/TSeries-1024cell-32concurrent-2reps-043.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-046.ty.h5"
tyh5Path = "/data/dlab/b115/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043.ty.h5"
fishDir = joinpath(splitpath(tyh5Path)[1:end-1]...)
expName = replace(splitpath(tyh5Path)[end], ".ty.h5" => "")
tseriesDir = joinpath(fishDir, expName)

tseries = h5read(tyh5Path, "/imaging/raw")
@assert size(tseries,4)==1
tseries = permutedims(tseries, (2,1,3,4,5))
tseries = tseries[:,:,:,1,:];

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

# read power
@assert length(glob("*.txt", fishDir)) == 1 # if not, need to be careful to choose
slmTxtFile = glob("*.txt", fishDir)[1]
stimGroupDF = CSV.File(open(read, slmTxtFile), header=["filepath", "powerFraction"]) |> DataFrame
stimGroupDF = stimGroupDF[trialOrder,:]

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
firstTargetGroup = matread.(findMatGroups(slmExpDir)[1])
lastTargetGroup = matread.(findMatGroups(slmExpDir)[end])
target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
# @assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}
is1024 = size(tseries,1)==1024
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
    is1024=is1024, zOffset=zOffset)

cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
@warn "hardcoded laser power"
slm1MaxPower = firstTargetGroup["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"]
slm1Power = 850mW
slm1powerPerCell = slm1Power * slm1MaxPower / 1000
cells[!, :slmNum] .= 1

cells[!, :laserPower] = round.(typeof(1.0mW), map(g->stimGroupDF.powerFraction[g], cells.stimGroup) .* slm1powerPerCell, digits=1)

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

winSize = Int(ceil(2*volRate))
# delay=Int(ceil(volRate*2))
delay=0

roiMask = constructROImasks(cells, H, W, Z, targetSizePx=targetSizePx);

cellsDF = makeCellsDF(tseries, cells, roiMask, winSize=winSize, delay=delay);
meanDf_f = combine(groupby(cellsDF, :cellID), :df_f => mean).df_f_mean
@show [sum(cellsDF.df_f[cellsDF.laserPower .== p] .> 0.4) for p in [3.6mW, 7.2mW, 10.8mW, 14.4mW]]
@show sum(meanDf_f .> 0.4)
cellIDrankings = sortperm(meanDf_f, lt=(>))
numNaN = sum(isnan.(cellsDF.df_f))
@assert numNaN == 0

@show sum(meanDf_f.>0.1), sum(meanDf_f.>0.4), sum(meanDf_f.>1.)

for power in sort(unique(cellsDF.laserPower))
    println("$power: $(mean(cellsDF.df_f[cellsDF.laserPower .== power]))")
end
# only 26 are >40% for two stim in control...
# @show first(cellsDF[rankings[numNaN+1:end],:],50)

##
# histogram(cellsDF.df_f)
Gadfly.plot(cellsDF, x="df_f",  ygroup="laserPower", Gadfly.Geom.subplot_grid(Gadfly.Geom.histogram))

##
before = Int(ceil(volRate*10))
after = Int(ceil(volRate*20))
# best
plots = []
nplots = 40

# nplots = 8
# for idx in rankings[end:-1:end-nplots+1]
# for idx in 1:nplots
for cellID in cellIDrankings[1:nplots]
# for idx in randperm(nplots)
    indices = findall(cellsDF.cellID .== cellID)
    push!(plots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after))
end
# p = plot(plots..., layout=(16, 4), legend=false, size=(1024,1024*2))
p = plot(plots..., layout=(8,5), size=(1024,1024*4))
# savefig("/home/tyler/Downloads/gcamp-control-random.png")
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