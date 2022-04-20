ENV["DISPLAY"] = "localhost:11"
##
using ImageView
using Lensman, PyCall, DataFrames, Gadfly, Distributed, StatsBase, Folds, ProgressMeter,
    Images, Statistics, HDF5, Base.Threads, LoopVectorization
using CSV
import PyPlot
import Lensman: @pun, @assign
import Unitful: ms
import Base.Threads: @threads
plt = PyPlot
matplotlib = plt.matplotlib
np = pyimport("numpy")
L = Lensman
##
resources = Resources()
recording = Recordings[
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
    # "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061"
    # "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
    # "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-118trial-122"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
    # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-lrhab-118trial-060"
    # "2021-06-15_rsChRmine_h2b6s_6dpf/fish1/TSeries-titration-192trial-050"
    # "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
](
    ; resources...,
    tseries_read_strategy=:lazy_tiff
);
# @pun tyh5_path = recording
##
@pun (tseries, voltageFile, stim_start_idx, stim_end_idx, frameStartIdx,
    tseriesH, tseriesZ, tseriesT, tseries_dir, frame_rate, fish_dir, exp_name,
    tseriesH, tseriesW, tseriesZ, tseriesT
) = recording;
##
frame_start_key="frame starts"
stim_key="respir"
sigma=3
voltages = CSV.File(open(read, voltageFile)) |> DataFrame
rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

# plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
frameStarted = diff(voltages[!,frame_start_key] .> std(voltages[!,frame_start_key]) .* sigma);
##
# add 1 for diff
frameStartIdx = findall(frameStarted .== 1) .+ 1
frame_duration = mean(diff(frameStartIdx))
ttlStarts = findTTLStarts(voltages[!,stim_key])
nstim_pulses = countPulsesMaxGap(ttlStarts)
stimStartIdx = ttlStarts[1:nstim_pulses:end]
ttlEnds = findTTLEnds(voltages[!,stim_key]);

##
# we want stim_start and stim_stop in (H,Z,T)
"Tested for 512x512 in B115."
function frame_frac_to_hzt(frame_frac, tseriesH, tseriesZ, scan_up=0.06)
    # add 1 for 1-indexing
    # frame 1 should be t=1; frame 10 should be t=1; frame 11 should be t=2
    T = Int(floor((frame_frac-1) / tseriesZ)) + 1
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
    d = Dict()
    @assign d = (startH,startZ,startT, endH,endZ,endT)
    return d
end
##

s = ttlStarts[rand(1:length(ttlStarts))]
res = get_ttl_hzt(s, frameStartIdx, tseriesH, tseriesZ, frame_rate)
@pun (startH,startZ,startT, endH,endZ,endT) = res


if (startH > 512) .& (endH > 512)
    println("both")
    ts = convert(Array{UInt16}, tseries[:,:,startZ,startT])
    im = imadjustintensity(ts)
    im = RGB.(im)
elseif startH > 512
    println("first")
    ts = convert(Array{UInt16}, tseries[:,:,endZ,endT])
    im = imadjustintensity(ts)
    im = RGB.(im)
    channelview(im)[1,1:endH,:] .= 0.3
elseif (endH > 512 ) | (startH > endH)
    println("second")
    ts = convert(Array{UInt16}, tseries[:,:,startZ,startT])
    im = imadjustintensity(ts)
    im = RGB.(im)
    channelview(im)[1,startH:tseriesH,:] .= 0.3
else
    println("neither")
    @show size(tseries[:,:,startZ,startT])
    ts = convert(Array{UInt16}, tseries[:,:,startZ,startT])
    im = imadjustintensity(ts)
    im = RGB.(im)
    channelview(im)[1,startH:minimum([endH,tseriesH]),:] .= 0.3
end
im
##
"For volume, remove horizontal band artifact that extends < 1 frame in height."
function remove_artifact_from_vol(vol, t, startH,startZ,startT,endH,endZ,endT; val=NaN)
    new_vol = copy(vol)
    if startT == t
        for z in 1:size(vol,3)
            if (z == startZ) && (z == endZ)
                # @info "artifact fully on this frame $z: $startH to $endH"
                new_vol[startH:minimum([endH,size(new_vol,1)]),:,z] .= val
            elseif (z == startZ)
                # @info "artifact starts on this frame"
                new_vol[startH:size(new_vol,1),:,z] .= val
            elseif (z == endZ)
                # @info "artifact ends on this frame"
                new_vol[1:endH,:,z] .= val
            end
        end
    elseif endT == t
        @assert endZ==1 "$endZ is not 1 for : $((startH,startZ,startT,endH,endZ,endT))"
        # @info "artifact ends on this frame"
        new_vol[1:endH,:,endZ] .= val
    end
    new_vol
end

function remove_artifacts_from_vol(vol, t, artifacts; val=NaN)
    idxs = (artifacts.startT .== t) .| (artifacts.endT .== t)
    relevant = artifacts[idxs,:]
    # new_vol =  convert(Array{Union{Missing,eltype(vol)}}, vol)
    if size(relevant,1) > 0
        # @info "removing $(size(relevant,1))"
        for row in eachrow(relevant)
            vol .= remove_artifact_from_vol(vol, t, row...; val=val)
        end
    end
    vol
end


function step_kalman(x,m; gain=0.8)
    # if ismissing(m)
    #     x
    # else
    #     x * gain + m * (1 - gain)
    # end
    isnan(m) ? x : x * gain + m * (1 - gain)
end

artifacts = map(ttlStarts) do s
    res = get_ttl_hzt(s, frameStartIdx, tseriesH, tseriesZ, frame_rate)
    @pun (startH,startZ,startT, endH,endZ,endT) = res
    startH,startZ,startT, endH,endZ,endT
end

cols = [:startH,:startZ,:startT, :endH,:endZ,:endT]
artifacts = DataFrame(artifacts)
rename!(artifacts, cols)

startT, startZ = artifacts[100,[:startT, :startZ]]
old_vol = tseries[:,:,:,startT]
new_vol = remove_artifacts_from_vol(old_vol, startT, artifacts; val=0)
new_vol = convert(Array{UInt16}, new_vol)
Gray.(imadjustintensity(new_vol[:,:,startZ]))

## Finally, save a denoised tseries!!

function kalman_filter_stream(tyh5_path, tseries, tseriesT, artifacts; save_uint16=true)
    h5 = h5open(tyh5_path, "w")
    tseries
    dset = create_dataset(h5, "kalman", datatype(UInt16), size(tseries))
    X = Float32.(tseries[:,:,:,1])
    @showprogress for t in 2:tseriesT
        vol = Float32.(tseries[:,:,:,t])
        M = remove_artifacts_from_vol(vol, t, artifacts)
        X = step_kalman.(X,M)
        if save_uint16
            # safe assuming max val of 8192; could overflow if 65000+
            dset[:,:,:,t] = convert(Array{UInt16}, round.(X))
        else
            dset[:,:,:,t] = X
        end
    end
    close(h5)
    tyh5_path
end
kalman_path = joinpath(fish_dir, exp_name * "_kalman.h5")
kalman_path = replace(kalman_path , "/scratch" => "/data/dlab")
mkpath(dirname(kalman_path))
# rm(kalman_path)
kalman_filter_stream(kalman_path, tseries, tseriesT, artifacts)



##
@pun tseries_dir = recording
h5 = h5open(kalman_path, "r")
dset = h5["kalman"]

# imshow(h5)

cmax = Int(round(quantile(tseries[:,:,4,60][:],0.95)))
Gray.(imadjustintensity(tseries[:,:,4,60]))
cmax = Int(round(quantile(dset[:,:,4,60][:],0.95)))
Gray.(imadjustintensity(dset[:,:,4,60], (0,cmax)))
close(h5)