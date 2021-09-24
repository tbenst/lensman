tifffile = pyimport("tifffile")

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
    if length(tiff_files) == 0
        error("No tif files found in $tifdir")
    end
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
    nstim_pulses = countPulsesMaxGap(ttlStarts)
    stimStartIdx = ttlStarts[1:nstim_pulses:end]
    stimEndIdx = findTTLEnds(voltages[!,stim_key])[nstim_pulses:nstim_pulses:end]
    stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimStartIdx]
    stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimEndIdx]
    stimStartFrameIdx, stimEndFrameIdx, frameStartIdx, nstim_pulses, ttlStarts
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
    idxs = parse.(Float64, map(x->x.attr["Index"], pvgp))
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
    cartIdxs = map(CartesianIndex ∘ Tuple, zip(y,x))
    # updated 2021-06-16 to return dict instead of array
    # no longer assumes a particular ordering of points
    Dict(Int(idx) => cartIdx for (idx,cartIdx) in zip(idxs, cartIdxs))
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
    @show names, idxs
    idxs_per_group = Dict(name => Int.(is) for (name,is) in zip(names,idxs))
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


# "Allows lazy but memoized indexing from HDF5 file."
# function lazy_read_tyh5(tyh5_path,
#         dset_path="/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm")
#     @info "assume WHZCT, convert to HWZT"
#     h5 = h5open(tyh5_path,"r", swmr=true)
#     dset = h5[dset_path]
#     return h5, LazyTy5(dset)
# end

function read_tyh5(tyh5_path,
        dset="/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm")
    @info "assume WHZCT, convert to HWZT"
    h5 = h5open(tyh5_path,"r")
    # dset = h5["/imaging/per_pixel_lstm_denoised"]
    # dset = h5["/imaging/per_pixel_lstm_denoised_maybe_longer_time"]
    dset = h5[dset]
    println(ndims(dset))
    if ndims(dset) == 5
        @assert size(dset)[4]==1
        W,H,Z,C,T = size(dset)
        dtype = eltype(dset)
        # if dtype == Float32
        #     # convert to UInt16, asume max value of 1...
        #     f_convert = float2uint
        # else
        #     f_convert = x -> x
        # end
        f_convert = x -> x
        tseries = zeros(UInt16, H, W, Z, T)
        # chunk by z for memory efficiency..
        @showprogress for z=1:Z
            # drop singleton channel
            dat = f_convert(dset[:,:,z,1,:])
            dat = permutedims(dat, (2,1,3))
            tseries[:,:,z,:] = dat
        end
    else
        error("stop")
        tseries = h5read(tyh5Path, "/imaging/per_pixel_lstm_denoised")
        # tseries = h5read(tyh5Path, "/imaging/raw")
        # @show size(tseries)
        # @assert size(tseries,2)==1
        # tseries = tseries[:,1,:,:,:]
        tseries = permutedims(tseries, (2,1,3,4))
        tseriesDir = joinpath(fishDir, expName)

        # tseries = h5read(tyh5Path, "/imaging/raw")
        # drop singleton channel
        # @assert size(tseries,4)==1
        # tseries = permutedims(tseries, (2,1,3,4,5))
        # tseries = tseries[:,:,:,1,:];
        # tseriesDir = joinpath(fishDir, expName)
    end
    tseries
end

"Read XYZ units in microns from Olympus .oir imaging file."
function read_oir_units(oir_file)
    x_regex = Regex(raw"^\s*<commonparam:x>([\d\.]+)+<")
    y_regex = Regex(raw"^\s*<commonparam:y>([\d\.]+)+<")
    z_regex = Regex(raw"^\s*<commonparam:z>([\d\.]+)+<")
    x_um = NaN
    y_um = NaN
    z_um = NaN
    for f in eachline(oir_file)
        x_match = match(x_regex, f)
        y_match = match(y_regex, f)
        z_match = match(z_regex, f)
        if ~isnothing(x_match)
            x_um = parse(Float64, x_match.captures[1])
        elseif ~isnothing(y_match)
            y_um = parse(Float64, y_match.captures[1])
        elseif ~isnothing(z_match)
            z_um = parse(Float64, z_match.captures[1])
        end
        if ~isnan(x_um) & ~isnan(y_um) & ~isnan(z_um)
            break
        end
    end
    (x_um, y_um, z_um)
end

"Glob one and only one file."
function glob_one_file(pattern, dir; nofail=false)
    files = glob(pattern, dir)
    if nofail
        if length(files) == 1
            return files[1]
        else
            return nothing
        end
    else
        @assert length(files)==1 "Expecting exactly one file for $dir/$pattern: $files"
        return files[1]
    end

end

function read_xml(xml_file)
    open(xml_file, "r") do io
        xml = read(io, String)
        xp_parse(xml)
    end
end

function load_zseries(zseries_dir; ch_str = "Ch3")
    tiff_files = joinpath.(zseries_dir,
        filter(x->(x[end-6:end]=="ome.tif") & occursin(ch_str, x),
        readdir(zseries_dir)))
    tif0 = ImageMagick.load(tiff_files[1])
    zseriesH, zseriesW = size(tif0)
    zseries_lateral_unit = microscope_lateral_unit(zseriesW)
    @info "hardcoded zseries z unit: 2μ"
    zseries_microscope_units = (zseries_lateral_unit, zseries_lateral_unit, 2.0μm)


    zseries = zeros(Normed{UInt16,16}, zseriesH, zseriesW, size(tiff_files, 1))
    @threads for z in 1:size(tiff_files,1)
        zseries[:,:,z] = ImageMagick.load(tiff_files[z])
    end
    AxisArray(zseries, (:y, :x, :z), zseries_microscope_units)
end

function read_suite2p(suite2p_dir)
    iscell = npzread(joinpath(suite2p_dir, "combined", "iscell.npy"));
    nCells = size(iscell,1)
    cells_mask = DataFrame()
    combined_stat = PyObject[]
    plane_dirs = glob("plane*", suite2p_dir)
    tseriesZ = length(plane_dirs)
    @assert tseriesZ > 0
    for z in 1:tseriesZ
        # have to process suite2p plane by plane due to
        # https://github.com/MouseLand/suite2p/issues/655
        stat = np.load(joinpath(suite2p_dir, "plane$(z-1)", "stat.npy"), allow_pickle=true);
        
        # iterate neurons in plane
        cell_idx_offset = length(combined_stat)
        for k in keys(stat)
            # uses 1 indexing (julia) rather than 0 indexing (suite2p/python)
            stat[k].update(iplane=z) # can't do []  assignment in PyCall
            # add ypix, xpix, cellNum to DataFrame
            xpix = stat[k].get("xpix")
            npix = length(xpix)
            cells_mask = vcat(cells_mask, DataFrame(xpix=xpix, ypix=stat[k].get("ypix"),
                              neuron_id=repeat([k+cell_idx_offset], npix)))
        end
        
        combined_stat = vcat(combined_stat, stat)
    end
    # nNeurons x 2
    cell_centers = aa2a(map(x->np.array(get(x,"med")), combined_stat));
    stat = Dict()
    @assign stat = (nCells, cell_centers, cells_mask, iscell)
    stat
end

function read_nwb_rois(nwb_path)
    ophys = h5open(nwb_path, "r")
    is_cell = transpose(read(ophys["/processing/ophys/ImageSegmentation/PlaneSegmentation/iscell"]))
    voxel_mask = read(ophys["/processing/ophys/ImageSegmentation/PlaneSegmentation/voxel_mask"])
    voxel_mask_index = read(ophys["/processing/ophys/ImageSegmentation/PlaneSegmentation/voxel_mask_index"])
    cell_traces = read(ophys["/processing/ophys/Fluorescence/Fluorescence/data"])
    cell_masks = DataFrame()
    for (i,e) in enumerate(voxel_mask_index)
        # https://hdmf.readthedocs.io/en/stable/hdmf.common.table.html#hdmf.common.table.VectorData
        if i == 1
            s = 1
        else
            s = voxel_mask_index[i-1] + 1
        end
        df = DataFrame(voxel_mask[s:e])
        df.cell_id = repeat([i], size(df,1))
        cell_masks = vcat(cell_masks, df)
    end
    nCells = size(is_cell,1)
    ret = Dict()
    @assign ret = (cell_traces, cell_masks, is_cell, nCells)
    ret
end

"Read HDF5 file in parallel and calculate per-trial df/f."
function get_df_f_per_voxel_per_trial_from_h5(tyh5_path, tseries_dset,
        stimStartIdx, stimEndIdx, window_len, tseriesH, tseriesW, tseriesZ)
    @info "rewrote to use LazyTy5; not tested"
    nTrials = length(stimStartIdx)
    df_f_per_voxel_per_trial = 
        SharedArray{Float64}(nTrials, tseriesH, tseriesW, tseriesZ);
    @sync @distributed for i in 1:nTrials
        proc_tseries = LazyTy5(tyh5_path, tseries_dset)
        s, e = stimStartIdx[i], stimEndIdx[i]
        f0 = mean(proc_tseries[:,:,:,s-window_len:s-1],dims=4)
        f = mean(proc_tseries[:,:,:,e+1:e+window_len],dims=4)
        df_f_per_voxel_per_trial[i,:,:,:] = @. (f - f0) / (f0 + 5)
    end
    convert(Array,df_f_per_voxel_per_trial)
end

function calc_save_trial_average(savepath::AbstractString, trial_average_dset::AbstractString,
    tseries, stim_start_idx, stim_end_idx, tseriesH, tseriesW, tseriesZ, trial_order;
    pre=window_len, post=window_len
)
    trial_average = calc_trial_average(tseries, stim_start_idx,
            stim_end_idx, tseriesH, tseriesW, tseriesZ, trial_order;
            pre=pre, post=post)
    h5write(savepath, trial_average_dset, trial_average)
    trial_average
end
    

"Calculate average response for each unique stimuli"
function calc_trial_average(tseries::LazyHDF5, stimStartIdx,
        stimEndIdx, tseriesH, tseriesW, tseriesZ, trialOrder;
        pre=16, post=16)
    nStimuli = maximum(trialOrder)
    nTrials = size(trialOrder, 1)
    nTrialsPerStimulus = [sum(trialOrder .== i) for i in 1:nStimuli]
    @assert nTrials == size(stimStartIdx, 1) # check TTL start-times match

    ## avg stim effec
    # TODO: should this be maximum instead to reduce chance of stim frame leaking in..?
    stim_trial_lens = stimEndIdx .- stimStartIdx
    @assert (maximum(stim_trial_lens) - minimum(stim_trial_lens)) <= 1
    trialTime = maximum(stim_trial_lens) + pre + post + 1
    
    # threads crashes HDF5
    stim_averages = map(1:nStimuli) do s
        stim_idxs = findall(trialOrder .== s)
        println("Calc average for sim $s")
        # segfaults..?
        # Folds.mapreduce((+), stim_idxs, DistributedEx()) do i
        @showprogress @distributed (+) for i in stim_idxs
            thread_tseries = typeof(tseries)(tseries.tyh5_path, tseries.dset_str)
            start = stimStartIdx[i]
            stop = start - pre + trialTime - 1
            trialType = trialOrder[i]
            arr = thread_tseries[:,:,:,start - pre:stop]
            close(thread_tseries)
            convert(Array{Float32}, arr)
        end
    end
    stim_averages =  stim_averages ./ nTrialsPerStimulus
    stim_averages = cat(stim_averages...; dims=5)
    permutedims(stim_averages, [1,2,3,5,4])
end

"Calculate average response for each unique stimuli"
function calc_trial_average(tseries::Union{Array{<:Real},LazyTiff}, stimStartIdx,
        stimEndIdx, tseriesH, tseriesW, tseriesZ, trialOrder;
        pre=16, post=16)
    nStimuli = maximum(trialOrder)
    nTrials = size(trialOrder, 1)
    nTrialsPerStimulus = [sum(trialOrder .== i) for i in 1:nStimuli]
    @assert nTrials == size(stimStartIdx, 1) # check TTL start-times match

    ## avg stim effec
    # TODO: should this be maximum instead to reduce chance of stim frame leaking in..?
    stim_trial_lens = stimEndIdx .- stimStartIdx
    @assert (maximum(stim_trial_lens) - minimum(stim_trial_lens)) <= 1
    trialTime = maximum(stim_trial_lens) + pre + post + 1
    # HWZCT

    
    stim_averages = map(1:nStimuli) do s
        stim_idxs = findall(trialOrder .== s)
        println("Calc average for stim $s")
        Folds.mapreduce((+), stim_idxs) do i
            start = stimStartIdx[i]
            stop = start - pre + trialTime - 1
            trialType = trialOrder[i]
            arr = tseries[:,:,:,start - pre:stop]
            convert(Array{Float32}, arr)
        end
    end
    stim_averages =  stim_averages ./ nTrialsPerStimulus
    stim_averages = cat(stim_averages...; dims=5)
    permutedims(stim_averages, [1,2,3,5,4])
end

"""
e.g.

read_zbrain_line("\$zbrain_dir/AnatomyLabelDatabase.hdf5",
    "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish")
"""
function read_zbrain_line(anatomy_label_h5_path, fishline, zbrain_units=zbrain_units;
        rostral=:left, dorsal=:down)
    volume = AxisArray(permutedims(
    h5read(anatomy_label_h5_path,
        fishline),
        (2,1,3)), (:y, :x, :z), zbrain_units)
    if rostral == :right
        volume = reverse(volume,dims=2);
    end
    if dorsal == :up
        volume = reverse(volume,dims=3);
    end
    volume
end

# "Read 4-channel olympus file."
# function read_olympus_btf(btf_file)
#     btf = collect(load(btf_file))
#     oZ = size(btf, 3)
#     # only blue for unknown reasons...
#     btf = reshape(btf,
#         size(btf)[1:2]..., 4, 2, Int(oZ//8))
#     # ch1 = 646nm (red)
#     channelview(btf)[1,:,:,1,:,:] = copy(channelview(btf)[3,:,:,1,:,:])
#     channelview(btf)[2:3,:,:,1,:,:] .= 0
#     # ch3 = 920nm (green)
#     channelview(btf)[2,:,:,3,:,:] = copy(channelview(btf)[3,:,:,3,:,:])
#     channelview(btf)[[1,3],:,:,3,:,:] .= 0
#     #ch4 = 405nm (blue), already okay
#     btf
# end

function read_olympus(oir_file)
    tifffile = pyimport("tifffile")
    x_um, y_um, z_um = read_oir_units(oir_file)
    @assert z_um == 1.0 # make sure not a eg 10um, quick-pan goof
    oir_img = tifffile.imread(replace(oir_file, "oir" => "ome.btf"));
    oir_img = permutedims(oir_img, [4,5,3,2,1]);
    AxisArray(oir_img, (:y, :x, :c, :t, :z),
        (y_um, x_um, 1, 1, z_um))
end

"Like joinpath, but make directory if doesn't exist."
function make_joinpath(parts...)
    path = joinpath(parts...)
    if ~isdir(path)
        mkpath(path)
    end
    path
end

"""
Using two channels to sync, write sparse datasets to HDF5.

Needed as HDF5 chokes with multithreading.
"""
function write_n_sparse_datasets(h5path, n, dset_size, channel,blocker)
    h5 = h5open(h5path, "w")
    h5["size"] = dset_size
    p = Progress(n; desc="Write masks:")
    while n > 0
        dset_name, data = take!(channel)
        take!(blocker)
        # save the sparse array
        # this used to sometimes error on all-zero arrays, but should
        # be fixed as of https://github.com/severinson/H5Sparse.jl/issues/5
        H5SparseMatrixCSC(h5, dset_name, data)
        n -= 1
        next!(p)
    end
    h5
end


function save_region_masks(region_mask_path, zseries, zbrain_masks, transform_paths,
    method=:cmtk;
    nprocs = 10, rostral=:right, dorsal=:up
)
    if method == :ants
        save_region_masks_ants(region_mask_path, zseries, zbrain_masks, transform_paths;
            nprocs = nprocs, rostral=rostral, dorsal=dorsal)
    elseif method == :cmtk
        save_region_masks_cmtk(region_mask_path, zseries, zbrain_masks, transform_paths;
            nprocs = nprocs, rostral=rostral, dorsal=dorsal)
    else
        error("unknown method $(string(method))")
    end
end

""""
Apply transform on Zbrain masks, and save sparse CSC arrays in h5 file.

ants_transforms: [affine, SyN] or [affine]

warping a single masks involves:
    1. read mask from database
    2. save mask as .nrrd
    3. apply the transform from ANTs (in my case) or CMTK (in your case) by calling the relevant shell command
    4. read in the warped result (likely a .nrrd file)
    5. save the warped result in a .h5 file

We want to do this for all masks--and there's ~300 of them. Fortunately, warping
using a known transform is pretty fast. Warping 300 masks with known transform is also
an embarrassingly parallel problem. Unfortunately, HDF5 as a format is a nightmare to
work on in parallel. the library segfaults unless extreme precautions are taken.
Basically, can only have one process write to a single .h5 file.

TO sovle, we use `nprocs` producers and a single consumer. Each producer writes the mask
name & the result of the warp into a channel. The consumer waits for the channel to fill,
and writes the result to the .h5 channel. Only the consumer has write access to the .h5
file.

The consumer is write_n_sparse_datasets and the producers are in the @distributed  block.
The bottleneck is actually the single consumer that writes to the .h5 file. this can
result in out of memory since 300 masks * 1024 * 1024 * 151 * sizeof(Float64) = 353 GB

Thus, a second channel is used, blocker. this channel has a capacity of `nprocs`,
which stops a consumer from proceeding if it is at capacity.

Approximately, we cap memory usage at nprocs * 1024 * 1024 * 151 * sizeof(Float64),
or say 21GB for 18 processes.
"""
function save_region_masks_ants(region_mask_path, zseries, zbrain_masks, ants_transforms;
    nprocs = 10, rostral=:right, dorsal=:up
)
    H, W, Z = size(zseries)
    mask_names = zbrain_masks["MaskDatabaseNames"]
    # nMasks = 10
    nMasks = length(mask_names)
    to_write = RemoteChannel(()->Channel(Inf))
    blocker = RemoteChannel(()->Channel(nprocs))
    @distributed for i in 1:nMasks
        @async begin
            put!(blocker, i)
            assignment = _region_mask_ants(i, zseries, zbrain_masks, ants_transforms, rostral, dorsal)
            put!(to_write, assignment)
            assignment = nothing
            GC.gc()
            next!(p)
        end
    end
    # TODO: perhaps we should use threads for producers, and one process for the consumer
    # should be an easy switch @distributed => @Threads
    # and remotecall_fetch for write_n_sparse_datasets
    # https://docs.julialang.org/en/v1/manual/distributed-computing/#Multi-processing-and-Distributed-Computing
    # currently, 1 hour 15 minutes for all masks.
    h5 = write_n_sparse_datasets(region_mask_path, nMasks, [H,W,Z], to_write, blocker)
    @everywhere GC.gc()
    return h5
end

function _region_mask_ants(i, zseries, zbrain_masks, ants_transforms, rostral, dorsal)
    H, W, Z = size(zseries)
    # get rid of `/` so not interpreted as a h5 path
    name = replace(zbrain_masks["MaskDatabaseNames"][i], "/" => "_")

    mask = read_mask(zbrain_masks, i; rostral=rostral, dorsal=dorsal)
    mask = AxisArray(Float32.(mask), AxisArrays.axes(mask))
    mask = antsApplyTransforms(zseries, mask,
        ants_transforms...) .> 0
    # TODO: or we could save as a .mat file for more standard format..?
    # https://github.com/JuliaIO/MAT.jl/blob/3ed629c05f7261e86c0dde0869d265e99a265efb/src/MAT_HDF5.jl
    mask = sparse(reshape(mask,length(mask),1))
    name => mask
end

function save_region_masks_cmtk(region_mask_path, zseries, zbrain_masks, cmtk_transform_path;
    nprocs=10, rostral=:right, dorsal=:up
)
    H, W, Z = size(zseries)
    zseries_path = joinpath(tmppath, ANTsRegistration.write_nrrd(zseries))
    mask_names = zbrain_masks["MaskDatabaseNames"]
    # nMasks = 2
    nMasks = length(mask_names)
    region_outline_path = replace(region_mask_path ,"region_masks.h5" =>"region_outlines.h5")
    
    n = nMasks
    # for type in [:mask, :outline]
    # for type in [:outline]
    # we use morphogradient instead of outline as cleaner 1px line
    for type in [:mask]
        if type == :mask
            h5 = h5open(region_mask_path, "w")
            N = nMasks
        elseif type == :outline
            h5 = h5open(region_outline_path, "w")
            @warn "missing 2 outlines for Left/Right Hemisphere"
            N = nMasks - 2
        end
        h5["size"] = [H, W, Z]

        @showprogress for i in 1:N
            name, data = _region_mask_cmtk(i, zseries_path, zbrain_masks, cmtk_transform_path,
                rostral, dorsal; type=type)
            if length(data.nzval) == 0
                # all zero
                h5[name] = 0
            else
                H5SparseMatrixCSC(h5, name, data)
            end
        end
    end
    nothing
end

function _region_mask_cmtk(i, zseries_path, zbrain_masks, cmtk_transform_path, rostral, dorsal;
    type=:mask
)
    # get rid of `/` so not interpreted as a h5 path
    name = replace(zbrain_masks["MaskDatabaseNames"][i], "/" => "_")

    mask = read_mask(zbrain_masks, i; rostral=rostral, dorsal=dorsal, type=type)
    mask = AxisArray(Float32.(mask), AxisArrays.axes(mask))

    # (first writing zseries to nrrd anyways so just pass path in)
    mask = apply_cmtk_transform(zseries_path, mask, cmtk_transform_path)
    mask = mask.data .> 0  # This thing visualizes correctly, so should be good

    # TODO: or we could save as a .mat file for more standard format..?
    # https://github.com/JuliaIO/MAT.jl/blob/3ed629c05f7261e86c0dde0869d265e99a265efb/src/MAT_HDF5.jl
    mask = sparse(reshape(mask, length(mask), 1))
    (name, mask)
end

function read_registered_mask(region_masks_h5, name; outline=false, verbose=false)
    if verbose
        @show name
    end
    shape = region_masks_h5["size"][:]
    ret = reshape(
        Array(sparse(H5SparseMatrixCSC(region_masks_h5, name))),
        # TODO benchmark
        # collect(sparse(H5SparseMatrixCSC(region_masks_h5, name))),
        shape...
    )
    if outline
        # only on XY dim
        return morphogradient(dilate(rm,[1,2]), [1,2])
    else
        return ret
    end        
end

"""Read all zbrain masks into sparse arrays.

Optionally read Tyler's hemisphere masks.
"""
function read_zbrain_masks(zbrain_dir; read_hemisphere=true)
    zbrain_masks = matread("$zbrain_dir/MaskDatabase.mat")
    zbrain_masks["MaskDatabaseNames"] = vec(zbrain_masks["MaskDatabaseNames"])
    if read_hemisphere
        left = sparse(H5SparseMatrixCSC("$zbrain_dir/hemisphere_masks.h5", "left"))
        right = sparse(H5SparseMatrixCSC("$zbrain_dir/hemisphere_masks.h5", "right"))
        push!(zbrain_masks["MaskDatabaseNames"], "Left hemisphere")
        push!(zbrain_masks["MaskDatabaseNames"], "Right hemisphere")
        zbrain_masks["MaskDatabase"] = hcat(zbrain_masks["MaskDatabase"], left, right)
        @warn "no outline for left/right hemisphere."
    end
    zbrain_masks["MaskDatabaseNames"] = replace.(zbrain_masks["MaskDatabaseNames"],
        "/" => "_")
    zbrain_masks
end

function read_cellpose_dir(directory)
    files = glob("MASK_plane_*.tif", directory)
    tifs = []
    startnum = 0
    for f in files
        tif = convert(Array{UInt32}, reinterpret(UInt16, ImageMagick.load(f)))
        nonzero_idxs = tif .> 0
        tif[nonzero_idxs] .+= startnum
        startnum = maximum(tif)
        push!(tifs,tif)
    end
    cat(tifs...; dims=3)
end
