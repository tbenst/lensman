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

function find_artifacts(ttlStarts, frameStartIdx, tseriesH, tseriesZ, frame_rate)
    artifacts = map(ttlStarts) do s
        res = get_ttl_hzt(s, frameStartIdx, tseriesH, tseriesZ, frame_rate)
        @pun (startH,startZ,startT, endH,endZ,endT) = res
        startH,startZ,startT, endH,endZ,endT
    end
    cols = [:startH,:startZ,:startT, :endH,:endZ,:endT]
    artifacts = DataFrame(artifacts)
    rename!(artifacts, cols)

end


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