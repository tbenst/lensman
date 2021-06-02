
"""Make .mat file for Sean's multiSLM stim software.

We minimize number of targetGroups for speed.
for Y shared: localToRemote = matpath -> "Y:\\" * replace(matpath[14:end], "/" => "\\")
"""
# ::Array{Float64,2}
function create_targets_mat(targets, outname::String; slmNum=1, frequency=5,
        localToRemote = matpath -> "O:\\" * replace(matpath[14:end], "/" => "\\"),
        targets_center = nothing)

    @assert frequency % 5 == 0 # sean SLM operates on 200ms clock
    @assert frequency <= 250 # max new stim every 4ms (on one SLM)
    @assert (slmNum == 1) | (slmNum == 2)

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
    @assert all(targets .<= 512)

    out_mat["targets"] = Array{Any}(nothing, 1, 1)
    out_mat["targets"][1] = copy(targets) # ensure not an Adjoint type
    out_mat["zrange"] = [minimum(targets[:,3]) maximum(targets[:,3])]
    # channel 1 (red) and 2 (green), used for GUI preview

    ntargets = size(targets, 1)
    
    # targetsGroup alternates between SLM1&2 (eg Float32(N,3) and zeros(0,0)) and each one is 2ms apart (4ms between first two SLM1 patterns)
    # maximum (200ms of time effectively)
    # @assert numGroups <= 99
    @warn "frequency is only accurate for 5Hz base clock"
    if frequency == 5
        numGroups = 4
    elseif frequency == 10
        # stim at 1 and 49
        numGroups = 50
    elseif frequency > 10
        numGroups = 98 # maximum for Sean's code
    end

    out_mat["targetsGroup"] = Array{Any}(nothing,1,numGroups)
    
    for i in 2:2:numGroups
        # SLM2
        out_mat["targetsGroup"][i] = zeros(0, 0)
    end
    for i in 1:2:numGroups
        # SLM1
        out_mat["targetsGroup"][i] = zeros(0, 0)
    end
    

    numStims = Integer(frequency/5)
    isi = 49 / numStims # max of 99; we ignore slot 99...
    stimIndices = Integer.(2 .* floor.(collect(0:numStims-1) * isi)) .+ slmNum

    for i in stimIndices
        out_mat["targetsGroup"][i] = copy(targets)
    end

    if ~isnothing(targets_center)
        out_mat["targetsCenter"] = Array{Any}(nothing,1,numGroups)
        for i in 1:numGroups
            out_mat["targetsCenter"][i] = targets_center
        end
    end

    # println("construct 10Hz pattern ")
    # out_mat["targetsGroup"][50] = copy(targets) # stim at +100ms

    matwrite(matpath, stim_struct)
    # strip /mnt/deissero... this is needed since Matlab on windows will read the file
    return localToRemote(matpath)
end

# ::Array{String,1}
"""Create trials txt file with matpath and power.

Defaults to 1 repetition per mat with power fraction 1. Can add multiple
by specifying powers
"""
function create_trials_txt(targets_mats, outname::String;
        outdir="/mnt/deissero/users/tyler/slm/masks/",
        repeatForPowers=[1])
    today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
    today_dir = mkpath(joinpath(outdir, today_str))
    txtpath = joinpath(today_dir, "$outname.txt")
    
    txt_io = open(txtpath, "w")
    for power in repeatForPowers
        for matpath in targets_mats
            println(txt_io, matpath * "\t$power")
        end
    end
    close(txt_io)
    return txtpath
end

"""Create trials txt file with matpath & power per matpath.

Unlike other method, no repetitions (1 power per matpath).
"""
function create_trials_txt(targets_mats, powers, outname::String;
        outdir="/mnt/deissero/users/tyler/slm/masks/")
    today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
    today_dir = mkpath(joinpath(outdir, today_str))
    txtpath = joinpath(today_dir, "$outname.txt")
    
    txt_io = open(txtpath, "w")
    for (power, matpath) in zip(powers, targets_mats)
        println(txt_io, matpath * "\t$power")
    end
    close(txt_io)
    return txtpath
end


"""Outname is for local filesystem (build platform).

powers: fractional power for that iteration. length > 1 means repeat all
frequencies: must be multiple of 5Hz. length > 1 means repeat all

Perform cartesian product of powers and frequencies.

TODO: this assumes a 5Hz clock for sean's code...but sometimes this is set
to eg 10Hz instead..."""
function create_slm_stim(target_groups, outname::String; slmNum=1,
        localToRemote = matpath -> "O:\\" * replace(matpath[14:end], "/" => "\\"),
        powers=[1], frequencies=[5], targets_center_list=nothing)
    @assert length(powers) == length(frequencies)
    if isnothing(targets_center_list)
        # for easy enumeration
        targets_center_list = [nothing for i in 1:length(target_groups)]
    end
    # create mat files; save windows path
    groupPowers = []
    targetsMats = []
    for (freq,power) in zip(frequencies, powers)
        @showprogress for (i, (targets, tc)) in enumerate(zip(target_groups, targets_center_list))
            name = "$(outname)_group_$(i)_freq_$(freq)"
            targetsMatPath = create_targets_mat(targets, name,
                slmNum=slmNum, localToRemote=localToRemote,
                frequency=freq, targets_center=tc)
            push!(targetsMats, targetsMatPath)
            push!(groupPowers,power)
        end
    end
    trials_txt = create_trials_txt(targetsMats, groupPowers, outname)
    return targetsMats, trials_txt
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
    tiff_files = joinpath.(tifdir, filter(x -> (length(x) >= 7) & (x[end - 6:end] == "ome.tif") &
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


"From a single voltageFile"
function getStimTimesFromVoltages(voltageFile, Z::Int;
        frame_start_key="frame starts", stim_key="respir", sigma=3)
    voltages = CSV.File(open(read, voltageFile)) |> DataFrame
    rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

    # plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
    frameStarted = diff(voltages[!,frame_start_key] .> std(voltages[!,frame_start_key]) .* sigma)
    frameStartIdx = findall(frameStarted .== 1) .+ 1

    ttlStarts = findTTLStarts(voltages[!,stim_key])
    stimDur = countPulsesMaxGap(ttlStarts)
    stimStartIdx = ttlStarts[1:stimDur:end]
    stimEndIdx = findTTLEnds(voltages[!,stim_key])[stimDur:stimDur:end]
    stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimStartIdx]
    stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimEndIdx]
    stimStartFrameIdx, stimEndFrameIdx, frameStartIdx
end

"""Read multiple voltageFiles and concatenate.
TperTrial is a bit of a evil hack.
Could potentially read from the CycleXXXX in voltageFile instead."""
function getStimTimesFromVoltages(voltageFiles::Array{String,1}, Z::Int, TperTrial::Int;
        frame_start_key="frame starts", stim_key="respir", sigma=3)
    start_idx = 0
    stimStartIdx = Int64[]
    stimEndIdx = Int64[]
    frameStartIdx = Int64[]
    nStimsPerTrial = 0
    for vf in voltageFiles
        tempStartIdx, tempEndIdx, tempFrameIdx = getStimTimesFromVoltages(vf, Z,
            frame_start_key=frame_start_key, stim_key=stim_key)
        tempStartIdx .+= start_idx
        tempEndIdx .+= start_idx
        stimStartIdx = vcat(stimStartIdx, tempStartIdx)
        stimEndIdx = vcat(stimEndIdx, tempEndIdx)
        frameStartIdx = vcat(frameStartIdx, tempFrameIdx)
        start_idx += TperTrial
        nStimsPerTrial = max(length(tempStartIdx), nStimsPerTrial)
    end
    stimStartIdx, stimEndIdx, frameStartIdx
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

"Recursively print keys."
function printMatKeys(mat; level=0, max_level=10)
    if level > max_level
        return nothing
    end
    for key in keys(mat)
        if typeof(key) == String
            println("  "^level * "$key")
            printMatKeys(mat[key], level=level+1, max_level=max_level)
        end
    end
end

function parseXML(xmlPath)
    local xml
    open(xmlPath, "r") do io
        xml = read(io, String)
        xml = xp_parse(xml);
    end
    xml
end

function getImagingPockels(prairieViewXML)
    parse(Float64,
        prairieViewXML[xpath"/PVScan/PVStateShard[1]/PVStateValue[@key='laserPower']/IndexedValue[1]/@value"][1])
end

"Add to xml doc for Bruker sequential stim"
function createMarkPointElement(xSeriesRoot, point::Int; numSpirals=10,
        initialDelay=0.15, interPointDelay=98)
    @assert initialDelay >= 0.15
    mpe = new_child(xSeriesRoot, "PVMarkPointElement")
    set_attribute(mpe, "Repetitions", "$numSpirals") # this seems to be ignored - gpl file overrides
    set_attribute(mpe, "UncagingLaser", "Uncaging")
    set_attribute(mpe, "UncagingLaserPower", "0") # gpl specifies...
    set_attribute(mpe, "TriggerFrequency", "None")
    set_attribute(mpe, "TriggerSelection", "None")
    set_attribute(mpe, "TriggerCount", "1")
    set_attribute(mpe, "AsyncSyncFrequency", "None")
    set_attribute(mpe, "VoltageOutputCategoryName", "None")
    set_attribute(mpe, "VoltageRecCategoryName", "Current")
    set_attribute(mpe, "parameterSet", "CurrentSettings")
    
    # nested under <PVMarkPointElement>!
    gpe = new_child(mpe, "PVGalvoPointElement")
    set_attribute(gpe, "InitialDelay", "$initialDelay")
    set_attribute(gpe, "InterPointDelay", "$interPointDelay")
    set_attribute(gpe, "Duration", "170") # from exported file but ????
    set_attribute(gpe, "SpiralRevolutions", "$numSpirals")
    set_attribute(gpe, "AllPointsAtOnce", "False")
    set_attribute(gpe, "Points", "Point $point")
    set_attribute(gpe, "Indices", "$point" )
end

function regex_glob(pattern, directory)
    file_names = readdir(directory)
    mask = occursin.(pattern, file_names)
    joinpath.(directory, file_names[mask])
end

""""Read Prairie View markpoints gpl list.

To measure, create markpoints at (1,1) and (1024,1024), export to .gpl, and view.

empirical measurements (B115)
zoom    X               Y
1       4.17170481      4.636140912
1.3     3.209      3.566
1.5     2.78113654      3.090760608
"""
function read_gpl(gpl_path::String; args...)
    gpl_xml = open(gpl_path, "r") do io
        gpl_xml = read(io, String)
        xp_parse(gpl_xml)
    end
    read_gpl(gpl_xml; args...)
end

function read_gpl(gpl_xml::ETree; zoom=1., width=1024, height=1024,
                  maxX=4.17170481, maxY=4.636140912,room="B115")
    pvgp = gpl_xml[xpath"/PVGalvoPointList/PVGalvoPoint"]
    x = parse.(Float64, map(x->x.attr["X"], pvgp))
    y = parse.(Float64, map(x->x.attr["Y"], pvgp))
    maxX /= zoom
    maxY /= zoom
    x .+= maxX
    x ./= 2*maxX
    x .*= width
    y .+= maxY
    y ./= 2*maxY
    y .*= height
    x = Int.(round.(x, digits=0))
    y = Int.(round.(y, digits=0))
    if room=="B115"
        # need to reverse for B115...
        x .= width .- x
        y .= height .- y # for b115 only...?
    end
    map(CartesianIndex ∘ Tuple, zip(y,x))
end

function read_markpoint_groups(gpl_path::String; args...)
    gpl_xml = open(gpl_path, "r") do io
        gpl_xml = read(io, String)
        xp_parse(gpl_xml)
    end
    read_markpoint_groups(gpl_xml; args...)
end

function read_markpoint_groups(gpl_xml::ETree; args...)
    points = read_gpl(gpl_xml; args...)
    groups = gpl_xml[xpath"/PVGalvoPointList/PVGalvoPointGroup"]
    aa_of_string = split.(map(x->x.attr["Indices"], groups), ",")
    names = map(x->x.attr["Name"], groups)
    idxs = map(y->parse.(Float64, y), aa_of_string)
    @info "Assume `Point x` = index + 1 in markpoints"
    idxs_per_group = Dict(name => Int.(is) .+ 1 for (name,is) in zip(names,idxs))
    Dict(name => map(i->points[i],gi) for (name,gi) in idxs_per_group)
end

""""Read Prairie View markpoints xml Series.
"""
function read_markpoints_series(xml_path; width=1024, height=1024)
    xml_file = open(xml_path, "r") do io
        xml_file = read(io, String)
        xp_parse(xml_file)
    end
    pvgp = xml_file[xpath"/PVMarkPointSeriesElements/PVMarkPointElement/PVGalvoPointElement/Point"]
    x = parse.(Float64, map(x->x.attr["X"], pvgp))
    y = parse.(Float64, map(x->x.attr["Y"], pvgp))
    x = Int.(round.(x*width, digits=0))
    y = Int.(round.(y*height, digits=0))
    map(CartesianIndex ∘ Tuple, zip(y,x))
end

getMatStimFreq(mat) = sum((~).(sum.(mat["cfg"]["exp"]["targets"]) .≈ 0.0))*5
getSLMnum(mat) = size(mat["cfg"]["exp"]["targets"][1]) == (0,0) ? 2 : 1




"""Write prairie view markpoints (galvo) gpl file.

neuron_locs is array of (x,y)."""
function write_markpoints(neuron_locs::Vector{CartesianIndex{2}}, filepath;
        W=512, magicX=7.6, magicY=8.3,spiral_size=0.228)
        # spiral_size=0.228 is 7um for B113
    # create an empty XML document
    xmlGpl = XMLDocument()
    # create & attach a root node
    xGplRoot = create_root(xmlGpl, "PVGalvoPointList")
    for (i,ci) in enumerate(neuron_locs)
        # x,y = Tuple(ci)
        y,x = Tuple(ci)
        # empirically, X ∈ [-7.6, 7.6], Y ∈ [8.3, -8.3]
        # where "(1,1)" is (-7.6,8.3)
        # and "(512,512)" is (7.6,-8.3)
        # and "(1,74)" meaning row 1 col 74 is (-5.37,8.3) 
        # and "(256,256)" is (0,0)
        # these magic numbers may differ by Bruker scope
        x = (x-W/2)/(W/2) # map to (-1,1)
        x *= magicX
        y = (y-W/2)/(W/2)
        y *= magicY
        
        gplList = new_child(xGplRoot, "PVGalvoPoint")
        # todo use Dict or named argument instead..?
        set_attribute(gplList, "X", "$x")
        set_attribute(gplList, "Y", "$y")
        set_attribute(gplList, "Name", "Point $i")
        set_attribute(gplList, "Index", "$(i-1)")
        set_attribute(gplList, "ActivityType", "MarkPoints")
        set_attribute(gplList, "UncagingLaser", "Uncaging")
        set_attribute(gplList, "UncagingLaserPower", "0.76")
        set_attribute(gplList, "Duration", "2")
        set_attribute(gplList, "IsSpiral", "True")
        set_attribute(gplList, "SpiralSize", "$spiral_size")
        set_attribute(gplList, "SpiralRevolutions", "5")
        set_attribute(gplList, "Z", "0")
    end

    save_file(xmlGpl, filepath)
end

function write_trial_order(trial_order, outname)
    trialOrderDF = DataFrame(copy(hcat(collect(1:length(trial_order)), trial_order)'))
    CSV.write(outname*"_trialOrder.txt", trialOrderDF, header=false, delim="\t")
    println("wrote $(outname*"_trialOrder.txt")")
    println("be sure to modify mSLM/SetupFiles/Experiments/<TODAY>/trialOrder.txt")
end