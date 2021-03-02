##
ENV["DISPLAY"] = "localhost:16.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow
import Gadfly
using Unitful: μm, m, s, mW

##
offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
zOffset = offset * 1e6

# tyh5Path = "/mnt/deissero/users/tyler/b115/2021-01-18_chrmine_kv2.1_h2b6s_6dpf/fish1_chrmine/TSeries-1024cell-32concurrent-2reps-043.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-046.ty.h5"
# tifDir = "/mnt/deissero/users/tyler/b115/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-076"
# tifDir = "/mnt/deissero/users/tyler/b115/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-10zplane-077"
tifDir = "/mnt/deissero/users/tyler/b115/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060"
# tifDir = "/mnt/deissero/users/tyler/b115/2021-02-02_f1_h33r_GC6f_6dpf/fish1_nochrmine/TSeries-1024cell-32concurrent-5power-044"
# tyh5Path = "/scratch/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043.ty.h5"
# tyh5Path = "/data/dlab/b115/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043.ty.h5"
# tyh5Path = "/scratch/2021-01-25_gcamp6f_6dpf/fish1/TSeries-1024cell-32concurrent-4power-044.ty.h5"

# TODO 1/29
# tyh5Path =  "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish1"
# tyh5Path =  "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish2/TSeries-1024cell-32concurrent-4power-043.ty.h5"
# tyh5Path =  "/mnt/b115_data/tyler/2021-01-25_rsChrmine_6f_6dpf/fish3/TSeries-1024cell-32concurrent-4power-046.ty.h5"
# tyh5Path =  "/data/dlab/b115/2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043.ty.h5"

# fishDir = joinpath(splitpath(tyh5Path)[1:end-1]...)
# expName = replace(splitpath(tyh5Path)[end], ".ty.h5" => "")

fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
expName = replace(splitpath(tifDir)[end], ".ty.h5" => "")
tylerSLMDir = joinpath(fishDir, "slm")


tseriesDir = joinpath(fishDir, expName)

# tseries = h5read(tyh5Path, "/imaging/raw")
# @assert size(tseries,4)==1
# tseries = permutedims(tseries, (2,1,3,4,5))
# tseries = tseries[:,:,:,1,:];

tseries = loadTseries(tifDir)
(H, W, Z, T) = size(tseries)

##
@show (H, W, Z, T)
slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

# TODO: wrap more of cell in function to reduce inline code..?
voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
2
# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
numETLvals = size(etlVals,1)
@assert numETLvals == Z
volRate = frameRate / Z
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

# read power
@assert length(glob("*.txt", tylerSLMDir)) == 1 # if not, need to be careful to choose

slmTxtFile = glob("*.txt", tylerSLMDir)[1]
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
##
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
##
winSize = Int(ceil(5*volRate))
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
open(tseriesDir*"_cellsDF.arrow", "w") do io
    tempDF = copy(cellsDF)
    tempDF.laserPower /= 1mW # can't serialize units...
    Arrow.write(io, tempDF)
end
## map of high df/f (top 50% of laser power only)
maxPower = maximum(cellsDF.laserPower)
Gadfly.with_theme(:dark) do
    mean_df = combine(groupby(cellsDF[cellsDF.laserPower .> mean(cellsDF.laserPower),:], :cellID), :df_f => mean)
    # mean_df = combine(groupby(cellsDF[cellsDF.laserPower .== maxPower,:], :cellID), :df_f => mean)
    mean_df = innerjoin(unique(cellsDF[:,[:x,:y,:cellID]]), mean_df, on=:cellID)
    dfmap = Gadfly.plot(mean_df, x=:x, y=:y, Gadfly.Geom.point, color=:df_f_mean,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true,aspect_ratio=1.))
    img = SVG(joinpath(plotDir,"$(expName)_dfmap_top_half_powers.svg"), 6inch, 5inch)
    Gadfly.draw(img, dfmap)
    dfmap
end


##

# histogram(cellsDF.df_f)
Gadfly.plot(cellsDF, x="df_f",  ygroup="laserPower", Gadfly.Geom.subplot_grid(Gadfly.Geom.histogram))

##
before = Int(ceil(volRate*10))
after = Int(ceil(volRate*20))
# best
plots = []
# nplots = 15
nplots = 45

# nplots = 8
# for idx in rankings[end:-1:end-nplots+1]
# for idx in 1:nplots
for cellID in cellIDrankings[1:nplots]
# for idx in randperm(nplots)
    indices = findall(cellsDF.cellID .== cellID)
    push!(plots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after))
end
# p = plot(plots..., layout=(16, 4), legend=false, size=(1024,1024*2))
# p = plot(plots..., layout=(5,3), size=(1000,800)) # 15 plots
p = plot(plots..., layout=(15,3), size=(1000,2400)) # 45 plots
savefig(p, joinpath(plotDir, "$(expName)_best45traces.png"))
# savefig(p, joinpath(plotDir, "$(expName)_best15traces_5sec_window.png"))
p
# savefig("/home/tyler/Downloads/gcamp-control-random.png")
##
# cellDFcontrol = Arrow.Table("/scratch/2021-01-25_gcamp6f_6dpf/fish1/TSeries-1024cell-32concurrent-4power-044_cellsDF.arrow") |> DataFrame
# arrowPath = "/mnt/deissero/users/tyler/b115/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-076_cellsDF.arrow"
arrowPath = tifDir*"_cellsDF.arrow"
df = Arrow.Table(arrowPath) |> DataFrame
df


##
# df = cellDFcontrol
df = copy(cellsDF)
df.laserPower /= 1mW


nreps = 1
n = maximum(df.cellID)*nreps
# idxs = sortperm(cellDFcontrol.df_f)
# Gadfly.plot(cellDFcontrol[idxs,:], x=:df_f, y=(1:n)./n, Geom.line, color=:laserPower)
p = nothing
# xmin = minimum(vcat(df.df_f, df.df_f))
# xmax = maximum(vcat(df.df_f, df.df_f))
xmin = -0.5
xmax = 2.0
cmax = maximum(df.laserPower)
for (i,pow) in enumerate(unique(df.laserPower))
    if i==1
        global p = plot(sort(df.df_f[df.laserPower .== pow],rev=true),
            (1:n)./n, line_z = pow, seriescolor=:lajolla, colorbar=false,
            label=pow, clims=(0,cmax), xlims=(xmin,xmax), ylabel="probability",
            xlabel="df/f")
    else
        p = plot!(sort(df.df_f[df.laserPower .== pow],rev=true),
            (1:n)./n, line_z = pow, seriescolor=:lajolla, colorbar=false,
            label=pow, clims=(0,cmax), xlims=(xmin,xmax))
    end
end
# savefig(p, joinpath(plotDir, "control_df_f_survival.png"))
savefig(p, joinpath(plotDir, "$(expName)_df_f_survival.png"))
p
##
n = length(cellDFcontrol.)
p = plot(sort(cellDFcontrol.df_f), (1:n)./n, 
    xlabel = "sample", ylabel = "Probability", 
    title = "Empirical Cumluative Distribution", label = "")

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