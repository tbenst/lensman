ENV["DISPLAY"] = "localhost:12"
##
using ImageView
using Lensman, PyCall, DataFrames, Gadfly, Distributed, StatsBase, Folds, ProgressMeter,
    Images, Statistics
using CSV
import PyPlot
import Lensman: @pun, @assign
import Unitful: ms
plt = PyPlot
matplotlib = plt.matplotlib
np = pyimport("numpy")
L = Lensman
##
resources = Resources()
recording = Recordings["2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"](
    ; resources...,
    tseries_dset=nothing, lazy_tiff=false
);
##
@pun (tseries, voltageFile, stim_start_idx, stim_end_idx, frameStartIdx,
    tseriesH, tseriesZ, tseriesT, tseries_dir, frame_rate
) = recording;
##
frame_start_key="frame starts"
stim_key="respir"
sigma=3
voltages = CSV.File(open(read, voltageFile)) |> DataFrame
rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

# plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
frameStarted = diff(voltages[!,frame_start_key] .> std(voltages[!,frame_start_key]) .* sigma)
##
# add 1 for diff
frameStartIdx = findall(frameStarted .== 1) .+ 1
frame_duration = mean(diff(frameStartIdx))
ttlStarts = findTTLStarts(voltages[!,stim_key])
nstim_pulses = countPulsesMaxGap(ttlStarts)
stimStartIdx = ttlStarts[1:nstim_pulses:end]
ttlEnds = findTTLEnds(voltages[!,stim_key])

##
# we want stim_start and stim_stop in (H,Z,T)
"Tested for 512x512 in B115."
function frame_frac_to_hzt(frame_frac, tseriesH, tseriesZ, scan_up=0.06)
    # add 1 for 1-indexing
    T = Int(floor(frame_frac / tseriesZ)) + 1
    # frame 1 should be z=1; frame 10 should be z=10; frame 11 should be z=1
    Z = Int(floor((frame_frac-1) % tseriesZ)) + 1
    adjH = tseriesH * (1+scan_up)
    H = Int(floor((frame_frac % 1) * adjH))
    (H,Z,T)
end
"""
Scan up is fraction of frame height spent scanning from bottom to top.
0.3ms buffer possibly okay; 0.1ms too tight; 0.5ms for reasonable safety
"""
function get_ttl_hzt(s, frameStartIdx, tseriesH, tseriesZ, frame_rate;
    scan_up=0.06, stim_duration=2ms, buffer=0.5ms
)
    frame_start = searchsortedlast(frameStartIdx, s)
    frame_duration = frameStartIdx[frame_start+1] - frameStartIdx[frame_start]
    frac_frame = Float64(s - 1 - frameStartIdx[frame_start]) / frame_duration
    frame_frac = frame_start + frac_frame
    
    frame_time = 1000ms/frame_rate
    before_frac = buffer / frame_time
    after_frac = (stim_duration + buffer) / frame_time


    startH, startZ, startT = frame_frac_to_hzt(frame_frac - before_frac,
        tseriesH, tseriesZ, scan_up)
    endH, endZ, endT = frame_frac_to_hzt(frame_frac + after_frac,
        tseriesH, tseriesZ, scan_up)

    return (startH,startZ,startT), (endH,endZ,endT)
end

s = ttlStarts[rand(1:length(ttlStarts))]
first, first_end = get_ttl_hzt(s, frameStartIdx, tseriesH, tseriesZ, frame_rate)
@show (first_pulse_H, first_pulse_Z, first_pulse_T) = first
@show (firste_pulse_H, firste_pulse_Z, firste_pulse_T) = first_end


if (first_pulse_H > 512) .& (firste_pulse_H > 512)
    println("both")
    im = imadjustintensity(tseries[:,:,first_pulse_Z,first_pulse_T])
    im = RGB.(im)
elseif first_pulse_H > 512
    println("first")
    im = imadjustintensity(tseries[:,:,firste_pulse_Z,firste_pulse_T])
    im = RGB.(im)
    channelview(im)[1,1:firste_pulse_H,:] .= 0.3
elseif (firste_pulse_H > 512 ) | (first_pulse_H > firste_pulse_H)
    println("second")
    im = imadjustintensity(tseries[:,:,first_pulse_Z,first_pulse_T])
    im = RGB.(im)
    channelview(im)[1,first_pulse_H:tseriesH,:] .= 0.3
else
    println("neither")
    im = imadjustintensity(tseries[:,:,first_pulse_Z,first_pulse_T])
    im = RGB.(im)
    channelview(im)[1,first_pulse_H:minimum([firste_pulse_H,tseriesH]),:] .= 0.3
end
im
##
test = tseries[:,:,5,100];

##
test[5,5] = missing