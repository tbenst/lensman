##
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, ProgressMeter, LightXML, LibExpat
using LightXML, Random, Arrow
import Gadfly
using Unitful: μm, m, s, mW

##
fishDir = "/data/dlab/b113/2021-04-27_wtChrmine_6f_6dpf/fish1"

## SAVE TRIAL ORDER / xml file for MarkPointSeries
nPoints = 21
outsideBrainNum = 20
nReps = 1
numStim = 10 # Hz since 1s stim... HARDCODE
dwellTime = 2 # ms
@warn "only setup for 1s stim, need to check if numStim is used correctly..."
numSpirals = 10
laserPowerPockels = 400
repetitions = 1
@warn "hardcode volRate"
volRate = 15


@assert nNeurons * repetitions <= 256 "approximately too big for Prairie View..."

xmlSeries = XMLDocument() # create & attach a root node
xSeriesRoot = create_root(xmlSeries, "PVSavedMarkPointSeriesElements")
set_attribute(xSeriesRoot, "Iterations", "1")
set_attribute(xSeriesRoot, "IterationDelay", delay)



# next, we make mark point series (xml) file
for r in 1:repetitions
    for i in randperm(nPoints)
        createMarkPointElement(xSeriesRoot,outsideBrainNum,numSpirals=1, initialDelay=4000)
        interPointDelay = (1000-numStim*dwellTime)/numSpirals
        createMarkPointElement(xSeriesRoot,i, numSpirals=numSpirals, interPointDelay=interPointDelay)
    end
end

save_file(xmlSeries, joinpath(fishDir, "$(nNeurons)cell_sequential.xml"))


## ANALYZE
tseries_dir = glob("TSeries*", fishDir)[end]
expName = splitpath(tseries_dir)[end]
tseries = loadTseries(tseries_dir)
(H, W, TperTrial, nRepeats) = size(tseries)
avg = mean(tseries, dims=4)[:,:,:,1];

plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

##
dataFolders = splitpath(tseries_dir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, _ = getExpData(xmlPath)
volRate = frameRate / Z

##
xml_path = glob("*MarkPoints.xml", tseries_dir)[1]
# x_magic, y_magic = markpoints_magic_numbers("B113")
# markpoints = read_gpl(xml_path, width=W, height=H, zoom=1, maxX=x_magic, maxY=y_magic)
markpoints = read_markpoints_series(xml_path, width=W, height=H)
targetSizePx = 8 # kinda a guess, should read in...
mean_img = RGB.(imadjustintensity(mean(avg,dims=3)[:,:,1]))
channelview(mean_img)[[1,3],:,:] .= 0
imshow(addTargetsToImage(mean_img, cartIdx2Array(markpoints),
    targetSize=targetSizePx))

##
roiMasks = constructROImasks(neuron_locs, H, W, targetSizePx);
@assert length(roiMasks) == nStimsPerTrial
##

## extract and save traces for easy plotting later...
fluor = DataFrame(time=Float64[], f=Float64[], cellID=UInt32[])

timeRange = (1:TperTrial) ./ volRate
for i in 1:nStimsPerTrial
    roiM = roiMasks[i]
    fluorescentTrace = extractTrace(avg, roiM)
    for (t,f) in zip(timeRange, fluorescentTrace)
        push!(fluor, (time=t, f=f, cellID=i))
    end
end

fluor_arrow_path = tseries_dir*"_cells_fluorescence.arrow"
if ~isfile(fluor_arrow_path)
    open(fluor_arrow_path, "w") do io
        Arrow.write(io, fluor)
    end
end
last(fluor,5)

##

function rawToDf_f(fluor::DataFrame; before=50)
    fluorDF = []
    # operate over each cell's trace...
    for df in groupby(fluor, [:cellID])
        kalmanFilt = imageJkalmanFilter(df.f)
        f0 = mean(kalmanFilt[1:before])
        df = copy(df)
        df_f = @. (kalmanFilt - f0) / (f0 + 10)
        insertcols!(df, size(df,2), :df_f => df_f)
        select!(df, Not(:f))
        push!(fluorDF, df)
    end
    fluorDF = vcat(fluorDF...)
end
fluorDF = rawToDf_f(fluor);
##

p = Gadfly.plot(fluorDF, x=:time, y=:df_f, ygroup=:cellID, Gadfly.Geom.subplot_grid(Gadfly.Geom.line))
img = SVG(joinpath(plotDir, "$(expName)_avg_trace.svg"), 10inch, 25inch)
Gadfly.draw(img, p)
p







##en
(H, W, T, nRepeats) = size(tseries)
tseries = reshape(tseries, H, W, 1, nRepeats*T);
TperTrial = T
(H, W, Z, T) = size(tseries)
# tseries = reshape(tseries, H, W, TperTrial, nRepeats); # undo

@show (H, W, Z, T)
## 

voltageFiles = glob("*VoltageRecording*.csv", tseries_dir)
vf = voltageFiles[1]

start_idx = 0
stimStartIdx = Int64[]
stimEndIdx = Int64[]
frameStartIdx = Int64[]
nStimsPerTrial = 0
for vf in voltageFiles
    tempStartIdx, tempEndIdx, tempFrameIdx = getStimTimesFromVoltages(vf, Z,
        frame_start_key="Frame begin", stim_key="Photodiode signal")
    tempStartIdx .+= start_idx
    tempEndIdx .+= start_idx
    stimStartIdx = vcat(stimStartIdx, tempStartIdx)
    stimEndIdx = vcat(stimEndIdx, tempEndIdx)
    frameStartIdx = vcat(frameStartIdx, tempFrameIdx)
    start_idx += TperTrial
    nStimsPerTrial = max(length(tempStartIdx), nStimsPerTrial)
end
    
@show length(stimStartIdx)
@assert length(frameStartIdx)==T # fails!
##
# nseconds = 3
# pre = Int(ceil(nseconds*volRate))+1
# post = Int(ceil(nseconds*volRate))+1

# trialOrder = repeat(collect(1:nStimsPerTrial), nRepeats)
# avgStim = trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder;
#     pre=pre, post=post);


##