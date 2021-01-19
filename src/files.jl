
"""Make .mat file for Sean's multiSLM stim software.

We minimize number of targetGroups for speed"""
# ::Array{Float64,2}
function create_targets_mat(targets, outname::String; slmNum=1)
    # outdir = "/mnt/deissero/users/tyler/slm/masks/")
    # today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
    # today_dir = mkpath(joinpath(outdir, today_str))
    # matpath = joinpath(today_dir, "$outname.mat")
    matpath = "$outname.mat"

    stim_struct = OrderedDict{String,Any}()
    # field order: im, targets, zrange, targetsGroup
    exp_key = "experimentS"
    # exp_key = "etb"
    stim_struct[exp_key] = Dict{String,Any}()
    out_mat = stim_struct[exp_key]

    out_mat["im"] = Array{Any}(nothing, 1, 2)
    # H, W = 1024, 1024
    H, W = 512, 512
    out_mat["im"][1] = zeros(H, W)
    out_mat["im"][2] = zeros(H, W)
    @assert all(targets .< 512)

    out_mat["targets"] = Array{Any}(nothing, 1, 1)
    out_mat["targets"][1] = copy(targets) # ensure not an Adjoint type
    out_mat["zrange"] = [minimum(targets[:,3]) maximum(targets[:,3])]
    # channel 1 (red) and 2 (green), used for GUI preview

    ntargets = size(targets, 1)
    
    # targetsGroup alternates between SLM1&2 (eg Float32(N,3) and zeros(0,0)) and each one is 2ms apart (4ms between first two SLM1 patterns)
    # maximum (200ms of time effectively)
    # @assert numGroups <= 99
    numGroups = 4
    out_mat["targetsGroup"] = Array{Any}(nothing,1,numGroups)
    
    for i in 2:2:numGroups
        # SLM2 is not used
        out_mat["targetsGroup"][i] = zeros(0, 0)
    end
    for i in 1:2:numGroups
        # future stim progressions go here
        out_mat["targetsGroup"][i] = zeros(0, 0)
    end
    
    if slmNum==1
        out_mat["targetsGroup"][1] = copy(targets) # stim at +0ms
    elseif slmNum==2
        out_mat["targetsGroup"][2] = copy(targets) # stim at +2ms
    else
        @error("bad slm num: $slmNum")
    end

    # println("construct 10Hz pattern ")
    # out_mat["targetsGroup"][50] = copy(targets) # stim at +100ms

    matwrite(matpath, stim_struct)
    # strip /mnt/deissero... this is needed since Matlab on windows will read the file
    return "O:\\" * replace(matpath[14:end], "/" => "\\")
end

# ::Array{String,1}
function create_trials_txt(targets_mats, outname::String;
        outdir="/mnt/deissero/users/tyler/slm/masks/")
    today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
    today_dir = mkpath(joinpath(outdir, today_str))
    txtpath = joinpath(today_dir, "$outname.txt")
    
    txt_io = open(txtpath, "w")
    for matpath in targets_mats
        println(txt_io, matpath * "\t1")
    end
    close(txt_io)
    return txtpath
end

"Outname is for local filesystem (build platform)"
function create_slm_stim(target_groups, outname::String; slmNum=1)
    targets_mats = []
    for (i, targets) in enumerate(target_groups)
        name = "$(outname)_group_$i"
        push!(targets_mats,create_targets_mat(targets, name, slmNum=slmNum))
    end
    trials_txt = create_trials_txt(targets_mats, outname)
    return targets_mats, trials_txt
end

# "Single-plane stim using sequential spirals."
# function create_sequential_bruker_stim(targetLocs2D; )
# end


function readZseriesTiffDir(tifdir; contains="Ch3")
    tiff_files = joinpath.(tifdir, filter(x -> (x[end - 6:end] == "ome.tif") & occursin(contains, x),
        readdir(tifdir)))
    tif0 = ImageMagick.load(tiff_files[1])
    H, W = size(tif0)

    zseries = zeros(Normed{UInt16,16}, H, W, size(tiff_files, 1))
    @threads for z in 1:size(tiff_files, 1)
        zseries[:,:,z] = ImageMagick.load(tiff_files[z])
    end
    zseries
end

"From bruker tif path, return (frameNum, planeNum)"
function getFramePlane(tifPath)
    m = match(r".*Cycle(\d+)_Ch\d_(\d+).ome.tif", tifPath)
    frameNum, planeNum = parse(Int, m[1]), parse(Int, m[2])
end


"""Prepare to read data in parallel.

Usage:
H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(".")
addprocs(18)
@everywhere using Lensman
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> _readTseriesTiffDir(S, framePlane2tiffPath, Z, T))

    """
function tseriesTiffDirMetadata(tifdir, containsStr="Ch3")
    tiff_files = joinpath.(tifdir, filter(x -> (x[end - 6:end] == "ome.tif") &
        occursin(containsStr, x), readdir(tifdir)))
    framePlane2tiffPath = Dict(getFramePlane(tp) => tp for tp in tiff_files)
    framePlanes = hcat(collect.(keys(framePlane2tiffPath))...)
    T = maximum(framePlanes[1,:])
    Z = maximum(framePlanes[2,:])
    if T == 1
        # single plane imaging, so T is last int in filename...
        framePlane2tiffPath = Dict((k[2], k[1]) => v for (k, v) in framePlane2tiffPath)
        T, Z = Z, T
    end
    tif0 = ImageMagick.load(tiff_files[1])
    H, W = size(tif0)
    H, W, Z, T, framePlane2tiffPath
end


""
function getStimTimesFromVoltages(voltageFile, Z::Int)
    voltages = CSV.File(open(read, voltageFile)) |> DataFrame
    rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

    # plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
    frameStarted = diff(voltages[!,"frame starts"] .> std(voltages[!,"frame starts"]) .* 3)
    frameStartIdx = findall(frameStarted .== 1) .+ 1

    ttlStarts = findTTLStarts(voltages[!,"respir"])
    stimDur = countPulsesMaxGap(ttlStarts)
    stimStartIdx = ttlStarts[1:stimDur:end]
    stimEndIdx = findTTLEnds(voltages[!,"respir"])[stimDur:stimDur:end]
    stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimStartIdx]
    stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimEndIdx]
    stimStartFrameIdx, stimEndFrameIdx
end


"""Parse Sean's filename to time
trialOrder_2020_11_02___21h03m29sr.txt -> Time(21,3,29)"""
function getTimeFromFilename(fn)
    re = r".*_(?<year>\d+)_(?<month>\d+)_(?<day>\d+)___(?<hour>\d+)h(?<min>\d+)m(?<sec>\d+)sr?.txt$"
    gs = match(re, fn)
    dt = parse.(Int, (gs["year"], gs["month"], gs["day"],
                     gs["hour"], gs["min"], gs["sec"]))
    DateTime(dt...)
end
