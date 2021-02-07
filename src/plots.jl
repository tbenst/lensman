using Plots, Plots.PlotMeasures
import Unitful: mW
function plotStim(tseries,roiMask,cells, idx::Int, volRate; before=30, after=60)
    roiM = roiMask[idx]
    x, y, z, _, stimStart, stimStop = cells[idx,:]
    theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
    theStart = maximum([stimStart-before, 1])
    plotRange = theStart:theEnd
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    fluorescentTrace = imageJkalmanFilter(fluorescentTrace)
    p = plot(plotRange/volRate, fluorescentTrace, left_margin=50px,
                ylims=(0,500))
        # xlabel="time (s)", ylabel="fluorescence")
    plot!(Lensman.rectangle((stimStop-stimStart)/volRate,maximum(fluorescentTrace),
        stimStart/volRate,0), opacity=.5, label="stim")
    p
end

"Plot multiple traces for the same cell on same plot."
function plotStim(tseries,roiMask,cells, indices::Array{Int,1}, volRate; before=30, after=60)
    idx = indices[1]
    roiM = roiMask[idx]
    stimStarts = cells[indices, :stimStart]
    stimStops = cells[indices, :stimStop]
    stimDur = maximum(stimStops .- stimStarts)
    # postTime = stimDur+after
    # timeRange = (-before+1:postTime)/volRate
    p = nothing
    maxF = 0.0
    cmax = maximum(cells.laserPower) / mW
    for idx in indices
        x, y, z, stimStart, stimStop, laserPower = cells[idx,[:x, :y, :z, :stimStart, :stimStop, :laserPower]]
        theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
        theStart = maximum([stimStart-before, 1])
        plotRange = theStart:theEnd
        timeRange = (plotRange .- stimStart) ./ volRate
        fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
        fluorescentTrace = imageJkalmanFilter(fluorescentTrace)
        if isnothing(p)
            p = plot(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=laserPower/mW, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
                # xlabel="time (s)", ylabel="fluorescence")
        else
            plot!(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=laserPower/mW, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
        end
        maxF = maximum([maxF, maximum(fluorescentTrace)])
    end
    plot!(Lensman.rectangle(stimDur/volRate,maxF,
        0,0), opacity=.5, label="stim")
    p
end

"Plot multiple traces for the same cell on same plot if recording multiple trials."
function plotStim(tseries,roiMask,cells, indices::Array{Int,1}, volRate, trial::Bool; before=30, after=60)
    if ~trial
        plotStim(tseries,roiMask,cells, indices, volRate; before=before, after=after)
    end
    idx = indices[1]
    roiM = roiMask[idx]
    stimStarts = cells[indices, :stimStart]
    stimStops = cells[indices, :stimStop]
    stimDur = maximum(stimStops .- stimStarts)
    # postTime = stimDur+after
    # timeRange = (-before+1:postTime)/volRate
    p = nothing
    maxF = 0.0
    cmax = maximum(cells.laserPower) / mW
    for idx in indices
        x, y, z, stimStart, stimStop, laserPower, trialNum = cells[idx,[:x, :y, :z, :stimStart, :stimStop, :laserPower, :trialNum]]
        theEnd = minimum([stimStop+after, size(tseries, ndims(tseries)-1)])
        theStart = maximum([stimStart-before, 1])
        plotRange = theStart:theEnd
        timeRange = (plotRange .- stimStart) ./ volRate
        fluorescentTrace = extractTrace(tseries[:,:,plotRange, trialNum], roiM)
        fluorescentTrace = imageJkalmanFilter(fluorescentTrace)
        if isnothing(p)
            p = plot(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=laserPower/mW, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
                # xlabel="time (s)", ylabel="fluorescence")
        else
            plot!(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=laserPower/mW, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
        end
        maxF = maximum([maxF, maximum(fluorescentTrace)])
    end
    plot!(Lensman.rectangle(stimDur/volRate,maxF,
        0,0), opacity=.5, label="stim")
    p
end