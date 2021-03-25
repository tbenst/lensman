using Plots, Plots.PlotMeasures
# import Gadfly
import Unitful: mW
function plotStim(tseries,roiMask,cells, idx::Int, volRate; before=30, after=60,
        colorBy=:laserPower)
    roiM = roiMask[cells[idx,:cellID]]
    x, y, z, _, stimStart, stimStop = cells[idx,:]
    theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
    theStart = maximum([stimStart-before, 1])
    plotRange = theStart:theEnd
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    fluorescentTrace = imageJkalmanFilter(medianfilt(fluorescentTrace))
    p = plot(plotRange/volRate, fluorescentTrace, left_margin=50px,
                ylims=(0,500))
        # xlabel="time (s)", ylabel="fluorescence")
    plot!(Lensman.rectangle((stimStop-stimStart)/volRate,maximum(fluorescentTrace),
        stimStart/volRate,0), opacity=.5, label="stim")
    p
end

"Plot multiple traces for the same cell on same plot."
function plotStim(tseries,roiMask,cells, indices::Array{Int,1}, volRate; before=30, after=60, colorBy=:laserPower)
    idx = indices[1]
    roiM = roiMask[cells[idx,:cellID]]
    stimStarts = cells[indices, :stimStart]
    stimStops = cells[indices, :stimStop]
    stimDur = maximum(stimStops .- stimStarts)
    # postTime = stimDur+after
    # timeRange = (-before+1:postTime)/volRate
    p = nothing
    maxF = 0.0
    cmax = maximum(cells[:,colorBy])
    if colorBy == :laserPower
        cmax = cmax / mW
    end
    for idx in indices
        x, y, z, stimStart, stimStop, color = cells[idx,[:x, :y, :z, :stimStart, :stimStop, colorBy]]
        theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
        theStart = maximum([stimStart-before, 1])
        plotRange = theStart:theEnd
        timeRange = (plotRange .- stimStart) ./ volRate
        fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
        fluorescentTrace = imageJkalmanFilter(medianfilt(fluorescentTrace))
        f0 = mean(fluorescentTrace[1:before])
        fluorescentTrace .-= f0
        fluorescentTrace ./= (f0 + 10)
        if colorBy == :laserPower
            color = color / mW
        end
        if isnothing(p)
            p = plot(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=color, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
                # xlabel="time (s)", ylabel="fluorescence")
        else
            plot!(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=color, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
        end
        maxF = maximum([maxF, maximum(fluorescentTrace)])
    end
    plot!(Lensman.rectangle(stimDur/volRate,maxF,
        0,0), opacity=.5, label="stim")
    p
end

"Plot multiple traces for the same cell on same plot if recording multiple trials."
function plotStim(tseries,roiMask,cells, indices::Array{Int,1}, volRate, trial::Bool; before=30, after=60, colorBy=:laserPower)
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
    cmax = maximum(cells[:, colorBy])
    if colorBy == :laserPower
        cmax = cmax / mW
    end
    for idx in indices
        x, y, z, stimStart, stimStop, color, trialNum = cells[idx,[:x, :y, :z, :stimStart, :stimStop, colorBy, :trialNum]]
        theEnd = minimum([stimStop+after, size(tseries, ndims(tseries)-1)])
        theStart = maximum([stimStart-before, 1])
        plotRange = theStart:theEnd
        timeRange = (plotRange .- stimStart) ./ volRate
        fluorescentTrace = extractTrace(tseries[:,:,plotRange, trialNum], roiM)
        fluorescentTrace = imageJkalmanFilter(medianfilt(fluorescentTrace))
        if colorBy == :laserPower
            color = color / mW
        end

        if isnothing(p)
            p = plot(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=color, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
                # xlabel="time (s)", ylabel="fluorescence")
        else
            plot!(timeRange, fluorescentTrace, left_margin=50px, legend=false,
                line_z=color, seriescolor=:lajolla, clims=(0,cmax), colorbar=false)
        end
        maxF = maximum([maxF, maximum(fluorescentTrace)])
    end
    plot!(Lensman.rectangle(stimDur/volRate,maxF,
        0,0), opacity=.5, label="stim")
    p
end

function max_95ci(x)
    confint(OneSampleTTest(x))[2]
end

function min_95ci(x)
    confint(OneSampleTTest(x))[1]
end

# chrmine_paper = Gadfly.Theme(
#     line_width=1mm,
#     minor_label_font="Arial",
#     major_label_font="Arial",
#     key_title_font="Arial",
#     key_label_font="Arial",
#     major_label_font_size=10pt,
#     minor_label_font_size=10pt,
#     key_title_font_size=10pt,
#     key_label_font_size=10pt
# )