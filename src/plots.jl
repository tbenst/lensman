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

mm2px = x -> Int(round(uconvert(u"inch",Quantity(x.value,u"mm"))*200/1u"inch",digits=0))

function inkscape_make_png(svg_fn, H, W;
        inkscape_path="/home/tyler/.nix-profile/bin/inkscape")
    png_fn = replace(svg_fn, ".svg" => ".png")
    px_w = mm2px(W)
    px_h = mm2px(H)
    
    cmd = `$inkscape_path -w $px_w -h $px_h $svg_fn --export-filename $png_fn`
    read(cmd, String)
end

"Symmetric log-like function."
function symlog(x,c=1/np.log(10))
    orig_sign = sign(x)
    out = log10(1 + abs(x/c))
    out * orig_sign
end

plus_sem(x) = mean(x)+sem(x)
minus_sem(x) = mean(x)-sem(x)
plus_std(x) = mean(x)+std(x)
minus_std(x) = mean(x)-std(x)

erodeN(n) = n == 1 ? erode : erode ∘ erodeN(n-1)
dilateN(n) = n == 1 ? dilate : dilate ∘ dilateN(n-1)
openingN(n) = dilateN(n) ∘ erodeN(n)

function opening_nonbool(x; N=1)
    nonzero = x .> 0
    nonzero = openingN(N)(nonzero)
    new = copy(x)
    new[(~).(nonzero)] .= 0
    new
end

opening_median = (x -> mapwindow(median!, x, (3,3))) ∘ opening_nonbool

rgb2tuple(rgb) = (red(rgb), green(rgb), blue(rgb))

"Remove salt & pepper."
function unsalt_n_pepper(im, felz_k=50, felz_min=10)
    im = opening_median(im)
    segments = felzenszwalb(im, felz_k, felz_min)
    im[segments.image_indexmap .== 1] .= 0
    im
end

"Remove salt & pepper."
function unsalt_n_pepper(vol::Array{T,3}; felz_k=50, felz_min=10) where T<:Real
    new = zeros(size(vol)...)
    for z in 1:size(new,3)
        im = opening_median(vol[:,:,z])
        segments = felzenszwalb(im, felz_k, felz_min)
        im[segments.image_indexmap .== 1] .= 0
        new[:,:,z] = im
    end
    new
end


mask2outline(m) = morphogradient(dilate(m,[1,2]), [1,2])