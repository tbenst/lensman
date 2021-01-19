using Plots, Plots.PlotMeasures
function plotStim(tseries,roiMask,cells, idx, volRate; before=30, after=60)
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