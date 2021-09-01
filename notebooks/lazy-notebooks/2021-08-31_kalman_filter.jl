ENV["DISPLAY"] = "localhost:12"
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
    # "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072"
](
    ; resources...,
    tseries_dset=nothing, lazy_tiff=true
);
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
@show (first_pulse_H, first_pulse_Z, first_pulse_T) = first
@show (firste_pulse_H, firste_pulse_Z, firste_pulse_T) = first_end


if (first_pulse_H > 512) .& (firste_pulse_H > 512)
    println("both")
    ts = convert(Array{UInt16}, tseries[:,:,first_pulse_Z,first_pulse_T])
    im = imadjustintensity(ts)
    im = RGB.(im)
elseif first_pulse_H > 512
    println("first")
    ts = convert(Array{UInt16}, tseries[:,:,firste_pulse_Z,firste_pulse_T])
    im = imadjustintensity(ts)
    im = RGB.(im)
    channelview(im)[1,1:firste_pulse_H,:] .= 0.3
elseif (firste_pulse_H > 512 ) | (first_pulse_H > firste_pulse_H)
    println("second")
    ts = convert(Array{UInt16}, tseries[:,:,first_pulse_Z,first_pulse_T])
    im = imadjustintensity(ts)
    im = RGB.(im)
    channelview(im)[1,first_pulse_H:tseriesH,:] .= 0.3
else
    println("neither")
    @show size(tseries[:,:,first_pulse_Z,first_pulse_T])
    ts = convert(Array{UInt16}, tseries[:,:,first_pulse_Z,first_pulse_T])
    im = imadjustintensity(ts)
    im = RGB.(im)
    channelview(im)[1,first_pulse_H:minimum([firste_pulse_H,tseriesH]),:] .= 0.3
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

##  NOT WORKING always
error("stop")
function vol_loader(channel,tseries, tseriesT, artifacts)
    println("vol_loader: init")
    @threads for t = 1:tseriesT
        println("vol_loader: $t")
        vol = Float32.(tseries[:,:,:,t])
        println("vol_loader: made vol")
        rem_vol = remove_artifacts_from_vol(vol, t, artifacts)
        println("vol_loader: blocked on $t")
        put!(channel, (t, rem_vol))
        println("vol_loader: sent $t")
    end
end

function vol_writer(tyh5_path, tseriesH, tseriesW, tseriesZ, tseriesT, denoised_vols)
    h5 = h5open(tyh5_path, "w")
    try
        dset = create_dataset(h5, "kalman", datatype(Float32), (tseriesH, tseriesW, tseriesZ, tseriesT))
        for t in 1:tseriesT
            println("vol_writer: blocked on $t")
            (tt,vol) = take!(denoised_vols)
            println("vol_writer: got $t")
            dset[:,:,:,tt] = vol
            # println("vol_writer: wrote $t")
        end
        close(h5)
    finally
        close(h5)
    end
    tyh5_path
end

function denoiser(raw_vols, denoised_vols, artifacts, tseriesT)
    println("denoiser: init")
    t, X = take!(raw_vols)
    first = Float32.(X)
    @assert t == 1
    println("denoiser: about to put")
    put!(denoised_vols, (1, first))
    println("denoiser: first put succeeded")
    # for t in 2:tseriesT
    for t in 950:tseriesT
        println("denoiser: blocked on take $t")
        tt, M = take!(raw_vols)
        println("denoiser: took $t")
        @assert t == tt
        # @threads for i in eachindex(X)
        # @sync for i in eachindex(X)
        #     @spawn  X[i] = step_kalman(X[i])
        # end
        # X = @tturbo step_kalman.(X,M)
        X = step_kalman.(X,M)
        println("denoiser: blocked on put $t")
        put!(denoised_vols, (t, convert(Array{Float32}, X)))
        println("denoiser: put $t")
    end
end

raw_vols = Channel(3 * Base.Threads.nthreads())
denoised_vols = Channel(3 * Base.Threads.nthreads())
# T = tseriesT
T = 61
loader_task = @spawn vol_loader(raw_vols, tseries, T, artifacts)
h5path = joinpath(fish_dir, exp_name * "_kalman.h5")
denoiser_task = @spawn denoiser(raw_vols, denoised_vols, artifacts, T)

kalman_path = vol_writer(joinpath(fish_dir, exp_name * "_kalman.h5"),
    tseriesH, tseriesW, tseriesZ, T, denoised_vols)











## Finally, save a denoised tseries!!

function kalman_filter_stream(tyh5_path, tseries, tseriesT, artifacts)
    h5 = h5open(tyh5_path, "w")
    tseries
    dset = create_dataset(h5, "kalman", datatype(Float32), size(tseries))
    X = Float32.(tseries[:,:,:,1])
    @showprogress for t in 2:tseriesT
        vol = Float32.(tseries[:,:,:,t])
        M = remove_artifacts_from_vol(vol, t, artifacts)
        X = step_kalman.(X,M)
        dset[:,:,:,t] = X
    end
    close(h5)
    tyh5_path
end
kalman_path = joinpath(fish_dir, exp_name * "_kalman.h5")
kalman_path = replace(kalman_path , "/scratch" => "/data/dlab")
kalman_filter_stream(kalman_path, tseries, tseriesT, artifacts)
# kalman_filter_stream(kalman_path, tseries, 6, artifacts)



##
h5 = h5open(kalman_path, "r")
dset = h5["kalman"]

# imshow(h5)

cmax = Int(round(quantile(tseries[:,:,4,60][:],0.95)))
Gray.(imadjustintensity(tseries[:,:,4,60]))
cmax = Int(round(quantile(dset[:,:,4,60][:],0.95)))
Gray.(imadjustintensity(dset[:,:,4,60], (0,cmax)))
close(h5)