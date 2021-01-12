module Lensman

using AxisArrays, ANTsRegistration, NIfTI, ImageMagick, Images,
    ImageDraw, ImageFiltering, PyCall, MAT, Dates, DataStructures,
    Statistics, SharedArrays, CSV, DataFrames, Suppressor, Plots,
    LinearAlgebra
import Base.Threads.@threads
using Distributed
import Unitful: μm
# https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules

include("algorithms.jl")

peak_local_max = PyNULL()
disk = PyNULL()
match_template = PyNULL()
np = PyNULL()
autolevel_percentile = PyNULL()
scipy = PyNULL()
machinery = PyNULL()
# py_utils = PyNULL()

function __init__()
    pyimport("skimage")
    copy!(peak_local_max, pyimport("skimage.feature").peak_local_max)
    copy!(disk, pyimport("skimage.morphology").disk)
    copy!(match_template, pyimport("skimage.feature").match_template)
    copy!(np, pyimport("numpy"))
    copy!(scipy, pyimport("scipy"))
    pyimport("scipy.ndimage")
    copy!(autolevel_percentile, pyimport("skimage.filters.rank").autolevel_percentile)

    # Nifty way to load local Python source code. This loads "src/py_utils" into a global `py_utils`
    copy!(machinery, pyimport("importlib.machinery"))
    # copy!(py_utils, pyimport("importlib.machinery").SourceFileLoader("py_utils", "src/py_utils.py").load_module("py_utils"))
end
    

include("Bruker.jl")
# include("segment.jl")

zbrain_units = (0.798μm, 0.798μm, 2μm)
tmppath = ANTsRegistration.userpath()

"Read (sparse) mask from Zbrain atlas and reshape to dense Array."
function read_mask(masks, idx; W=621, H=1406, Z=138, units=zbrain_units)
    example_mask = reshape(Array(masks["MaskDatabase"][:,idx]), (H, W, Z))
    # reshape as [621, 1406, 138]
    example_mask = permutedims(example_mask, (2, 1, 3))
    AxisArray(example_mask, (:y, :x, :z), units)
end

"Combine Zbrain masks whose name includes the given string."
function getMaskNameIncludes(masks, maskNameIncludes, units=zbrain_units)
    masks_raphe_idxs = findall(occursin.(maskNameIncludes, masks["MaskDatabaseNames"])) |>
        x -> map(ci -> ci[2], x)
    mask = foldl((a, n) -> a .| read_mask(masks, n), masks_raphe_idxs;
        init=read_mask(masks, masks_raphe_idxs[1]))
    mask = AxisArray(Float32.(mask), (:y, :x, :z), zbrain_units)
end

function antsApplyTransforms(fixedPath::String, moving::AxisArray, transformPath::String)
    movingPath = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
    antsApplyTransforms(fixedPath, movingPath, transformPath)
end

function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath::String)
    maskoutname = joinpath(tmppath, ANTsRegistration.randstring(10) * ".nii.gz")
    antsApplyTransforms(fixedPath, movingPath, transformPath, maskoutname)
end

function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath::String, maskoutname::String)
    run(`antsApplyTransforms --float -d 3 -i $movingPath -r $fixedPath -t $transformPath -o $maskoutname`)
    niread(maskoutname)
end

function read_zseries(tiff_files)
    tif0 = ImageMagick.load(tiff_files[1])
    H, W = size(tif0)
    zseries = zeros(Normed{UInt16,16}, H, W, size(tiff_files, 1))
    @threads for z in 1:size(tiff_files, 1)
        zseries[:,:,z] = ImageMagick.load(tiff_files[z])
    end
    return zseries
end

function create_nuclei_mask(microscope_units)
    nucleus_r = (2.75μm / microscope_units[1])μm
    space = 0.3μm # WARN: need to adjust CirclePointRadius, too
    cell_r = nucleus_r + space
    template_size = Int32(ceil(2 * (cell_r) / microscope_units[1]))
    # ensure odd size
    template_size += abs((template_size) % 2 - 1)
    # we make large and later resize for anti-aliasing
    nucleus_template = Gray.(zeros(513, 513))
    draw!(nucleus_template, Ellipse(CirclePointRadius(256, 256, 192)))
    nucleus_template = Float64.(nucleus_template)
    # nucleus_template[nucleus_template.==0.] .= -1.
    sz = (template_size, template_size)
    σ = map((o, n) -> 0.75 * o / n, size(nucleus_template), sz)
    kern = KernelFactors.gaussian(σ)

    nucleus_template = imresize(imfilter(nucleus_template, kern, NA()), sz)
end

"""Return center coordinates for each neuron.

TODO: perform in 3D?!
"""
function findNeurons(plane, thresh_adjust=1.2, featSize=2, maxiSize=4)
    # mask out background
    adj_plane = adjust_histogram(imadjustintensity(plane), GammaCorrection(0.5))
    threshold = otsu_threshold(adj_plane) * thresh_adjust
    otsu_mask = adj_plane .> threshold
    img = adj_plane .* otsu_mask
    feat_img = match_template(Float64.(img), disk(featSize), pad_input=true, mode="constant", constant_values=0)
    local_maxi = peak_local_max(feat_img, indices=false, footprint=disk(maxiSize), exclude_border=false, threshold_rel=.1)
end

"""Make .mat file for Sean's multiSLM stim software.

We minimize number of targetGroups for speed"""
# ::Array{Float64,2}
function create_targets_mat(targets, outname::String)
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
    
    # perhaps this is optional...?
    # targetsGroup alternates between SLM1&2 (eg Float32(N,3) and zeros(0,0)) and each one is 2ms apart (4ms between first two SLM1 patterns)
    # numGroups = 99 # not sure why this varies...
    numGroups = 4
    out_mat["targetsGroup"] = Array{Any}(nothing, 1, numGroups)
    out_mat["targetsGroup"][1] = copy(targets) # stim at +0ms

    for i in 2:2:numGroups
        # SLM2 is not used
        out_mat["targetsGroup"][i] = zeros(0, 0)
    end
    for i in 3:2:numGroups
        # future stim progressions go here
        out_mat["targetsGroup"][i] = zeros(0, 0)
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
function create_slm_stim(target_groups, outname::String)
    targets_mats = []
    for (i, targets) in enumerate(target_groups)
        name = "$(outname)_group_$i"
        push!(targets_mats, create_targets_mat(targets, name))
    end
    trials_txt = create_trials_txt(targets_mats, outname)
    return targets_mats, trials_txt
end


function cartIdx2Array(x::Vector{<:CartesianIndex})
    # make tall matrix
    permutedims(reduce(hcat, collect.(Tuple.(x))), (2, 1))
end

"ImageJ 'kalman filter'"
function imageJkalmanFilter(ys;gain=0.8,useBurnIn=true)
    predictions = zeros(size(ys))
    predictions[1] = ys[1]
    for t in 1:size(ys, 1) - 1
        K = useBurnIn ? 1 / (t + 1) : 0.
        pred = predictions[t]
        obs = ys[t + 1]
        predictions[t + 1] = pred * gain + obs * (1 - gain) + K * (obs - pred)
    end
    predictions
end

function extractTrace(tseries, mask::Array{Bool,T}) where T
    extractTrace(tseries, findall(mask))
end

function extractTrace(tseries, cartIdxs::Array{T,1}) where T <: CartesianIndex
    nT = size(tseries, ndims(tseries))
    trace = zeros(Float64, nT) 
    for t in 1:nT
        trace[t] = mean(reinterpret(UInt16, view(tseries, cartIdxs, t)))
    end
    trace
end

"Calculate df/f, assuming last dimension is time."
function df_f(f_timeseries, f0_timeseries)
    d = ndims(f_timeseries)
    f = selectdim(mean(f_timeseries, dims=d), d, 1)
    f0 = selectdim(mean(f0_timeseries, dims=d), d, 1)
    @. (f-f0)/f
end

cartesianIndToArray = cartIdx2Array

"Given (2d) image and list of CartesianIndices, return RGB image."
function colorImage(image, cartesianIdxs; nDilate=0, imAlpa=0.5)
    im = RGB.(image)
    stim_points = zeros(Bool, size(im))
    stim_points[cartesianIdxs] .= true
    for _ in 1:nDilate
        stim_points = dilate(stim_points)
    end
    channelview(im)[[1,3],:,:,:] .= 0
    channelview(im)[1,:,:,:] .= imAlpa .* float(stim_points)
    im
end


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


# function readTseriesTiffDir(tifdir::String; contains="Ch3")
#     tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") &
#         occursin(contains, x), readdir(tifdir)))
#     framePlane2tiffPath = Dict(getFramePlane(tp) => tp for tp in tiff_files)
#     framePlanes = hcat(collect.(keys(framePlane2tiffPath))...)
#     T = maximum(framePlanes[1,:])
#     Z = maximum(framePlanes[2,:])
#     tif0 = ImageMagick.load(tiff_files[1])
#     H, W = size(tif0)

#     # tseries = zeros(Normed{UInt16,16}, H, W, Z, T)
#     tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
#         init = S -> _readTseriesTiffDir(S, framePlane2tiffPath, Z, T))
# end

stripLeadingSpace(s) = (s[1] == ' ') ? s[2:end] : s

function findTTLStarts(voltages)
    frameStarted = diff(voltages .> std(voltages) .* 3)
    frameStartIdx = findall(frameStarted .== 1) .+ 1
end

function findTTLEnds(voltages)
    frameStarted = diff(voltages .< std(voltages) .* 3)
    frameStartIdx = findall(frameStarted .== 1) .+ 1
end


perm(permutation, array) = map(x -> permutation[x], array)

"""Pass image (HxW) and targetLocs (Nx2), and return image HxW"

targetSize=4.41 comes from 7um spiral on 25x objective @ 1x (0.63 μm/px)
"""
function addTargetsToImage(image::Array{RGB{T},2}, targetLocs::Array{<:Real,2};
  channel="red", targetSize::Real=11.11, alpha=0.4) where T <: Real
    stim = zeros(Gray, size(image))
    for (x, y) in eachrow(targetLocs)
        # draw!(stim, Ellipse(CirclePointRadius(y,x,targetSize)))
        draw!(stim, CirclePointRadius(y, x, targetSize))
    end
    
    if channel == "red"
        c = 1
    elseif channel == "blue"
        c = 3
    elseif channel == "red"
    c = 2
    else
        throw(ArgumentError("Bad channel $channel"))
    end
    channelview(image)[c,:,:,:] .= stim * alpha
    image
end


"Convert mask to an Array of indices, similar to argwhere in Python."
function mask2ind(mask)
    transpose(hcat(collect.(Tuple.(findall(mask)))...))
end

# /2: account for sean's code using 512 x 512 space
# while we image in 1024x1024
cartIdx2SeanTarget(ind, z) = hcat(reverse(collect(Tuple(ind)) / 2)..., z)

cartIdx2SeanTarget512(ind, z) = hcat(reverse(collect(Tuple(ind)))..., z)

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

"Count TTLs. TODO maxGap should be read from .xml..."
function countPulsesMaxGap(stimStartIdx, maxGap=1000)
    count = 1
    for i in 2:size(stimStartIdx, 1)
        if stimStartIdx[i] - stimStartIdx[i - 1] < maxGap
            count += 1
        else
            break
        end
    end
        count
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

function getStimFrameIdx(voltageFile, Z; colname="respir", nTTLperStim=10)
    @warn "deprecated"
    voltages = CSV.File(open(read, voltageFile)) |> DataFrame
    rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)]);

    # ! means get reference for df, not copy
    frameStarted = diff(voltages[!,"frame starts"] .> std(voltages[!,"frame starts"]) .* 2)
    frameStartIdx = findall(frameStarted .== 1) .+ 1;

    # eg for 10vol stim (2s @ 5Hz), nTTLperStim=10
    stimStartIdx = findTTLStarts(voltages[!,"respir"])[1:nTTLperStim:end - nTTLperStim + 1]
    stimEndIdx = findTTLEnds(voltages[!,"respir"])[nTTLperStim:nTTLperStim:end];
    ##
    stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimStartIdx]
    stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s) / Z)) for s in stimEndIdx]
    stimStartFrameIdx, stimEndFrameIdx
end

function getExpData(xmlPath)
    local tseries_xml
    open(xmlPath, "r") do io
        tseries_xml = read(io, String)
        tseries_xml = xp_parse(tseries_xml)
    end;
    # warning: find has different syntax from xpath!!!
    expDate = LibExpat.find(tseries_xml, "/PVScan[1]{date}")
    expDate = DateTime(expDate, "m/d/y I:M:S p")
    # much faster by prepending /PVScan/PVStateShard[1] to search first frame only
    framePeriod = parse(Float64, tseries_xml[xpath"/PVScan/PVStateShard[1]/PVStateValue[@key='framePeriod']/@value"][1])
    frameRate = 1 / framePeriod
    etlVals = getPlaneETLvals(tseries_xml)
    expDate, frameRate, etlVals
end

"""
only times that are within subfolders were used by a "Run"
EXAMPLE: slmExpDir=/mnt/deissero/users/tyler/b115/SLM_files/02-Nov-2020
e.g. after generating holograms, on first run (21h09m31s),
a subfolder is created like 2020_11_2_21_6_50_, and Sean's code 
copies trialOrder_2020_11_02___21h09m31s.txt to it.

If we "Stop" and "Run" again without regenerating (but possibly
changing trial order), no new folder is made BUT the new files for
the RUN are copied to the same subfolder
(trialOrder_2020_11_02___22h04m47s.txt)
"""
function getTrialOrder(slmExpDir, expDate)    
    slmRunDirs = joinpath.(slmExpDir,
        filter(x -> isdir(joinpath(slmExpDir, x)), readdir(slmExpDir)))

    trialOrderTxts = vcat(glob.("trialOrder_*s.txt", slmRunDirs)...)
    tempDts = getTimeFromFilename.(trialOrderTxts)

    # sort by datetime
    trialOrderTxts = trialOrderTxts[sortperm(tempDts)]
    availableDateTimes = getTimeFromFilename.(trialOrderTxts)

    # slmDate should always be younger than expDate
    dtIdx = searchsortedfirst(availableDateTimes, expDate) - 1
    slmDateTime = availableDateTimes[dtIdx]
    timeDiff = Dates.Second(expDate - slmDateTime)
    println("started SLM $(timeDiff) seconds before imaging")
    if (timeDiff > Second(120)) | (timeDiff < Second(0))
        @warn "timeDiff looks funky...verify we have correct trial"
    end
    trialOrderTxt = trialOrderTxts[dtIdx]

    trialOrder = CSV.File(trialOrderTxt, delim="\t", header=true) |> Tables.matrix
    @assert size(trialOrder, 1) == 1
    trialOrder[1,:], joinpath(splitpath(trialOrderTxt)[1:end - 1]...)
end

function loadTseries(tifdir)
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir)
    tseries = Array{UInt16}(undef, H, W, Z, T)
    p = Progress(T, 1, "Load Tseries: ")
    
    @threads for t in 1:T
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            tseries[:,:,z,t] .= reinterpret(UInt16, ImageMagick.load(tp))
        end
        next!(p)
    end
    tseries
end

function write_experiment_to_tyh5(fish_dir, exp_name, output_path)
    # Convention on output path?
    tif_dir = joinpath(fish_dir, exp_name);
    tseries = loadTseries(tif_dir);

    loader = machinery.SourceFileLoader("py_utils", "src/py_utils.py")
    py_utils = loader.load_module("py_utils")
    py_utils.write_tseries_to_tyh5(tseries, output_path)
    size(tseries)
    
    # TODO: write out stim pattern if it's available
    #   - tricky with the zOffset and such

    # TODO: Is this seriously not enough?
    tseries = nothing
    GC.gc()
end

function trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder; pre=16, post=16)
    nStimuli = maximum(trialOrder)
    nTrials = size(trialOrder, 1)
    nTrialsPerStimulus = Int(size(trialOrder, 1) / nStimuli)
    @assert nTrials == size(stimStartIdx, 1) # check TTL start-times match 

    ## avg stim effec
    trialTime = minimum(stimEndIdx .- stimStartIdx) + pre + post + 1
    # HWZCT
    avgStim = zeros(size(tseries)[1:3]..., nStimuli, trialTime);
    ## sloooow
    for (i, start) in enumerate(stimStartIdx)
        stop = start - pre + trialTime - 1
        trialType = trialOrder[i]
        avgStim[:,:,:,trialType,:] .+= tseries[:,:,:,start - pre:stop]
    end
    avgStim ./= nTrialsPerStimulus;
end

sym_adjust_gamma(image, gamma) = sign.(image) .* adjust_gamma(abs.(image), gamma)

function constructGroupMasks(groupLocs, H, W, Z;
  targetSizePx=(7μm * 14.4 / 25) / 0.6299544139175637μm)
    groupMasks = Gray.(zeros(Bool, H, W, Z, size(groupLocs, 1)))
    for (g, group) in enumerate(groupLocs)
        for (x, y, z) = eachrow(group)
            draw!(view(groupMasks, :, :, z, g), Ellipse(CirclePointRadius(x, y, targetSizePx)))
        end
    end
    reinterpret(Bool, groupMasks)
end


function init_workers(nprocs=36) 
    addprocs(nprocs)
    exp = quote
        @everywhere begin
            import Pkg
            Pkg.activate(".")
            using Lensman
        end
    end
    @suppress_err eval(macroexpand(Distributed, exp))
end

function allWithin(array, fraction)
    theMean = mean(array)
    threshold = ceil(fraction * theMean)
    @assert all(abs.(array .- theMean) .< threshold)
end

function findMatGroups(outdir;
        group_re=Regex("param_BHV004_(\\d+)\\.mat\$"))
    mats = glob("*.mat", outdir)
    # ensure we match the name/regex
    matMatches = filter(m -> m != nothing, match.(group_re, mats))
    mats = joinpath.(outdir, getfield.(matMatches, :match))
    group_nums = parse.(Int, getindex.(matMatches, 1))
    mats = mats[sortperm(group_nums)]
end

function rpadArray(array, n)
    toPad = n - size(array, 1)
    parent(padarray(array, Fill(0, (0,), (toPad,))))
end

"Create DataFrame with row per individual cell stimulated."
function makeCellsDF(target_groups, stimStartIdx, stimEndIdx, trialOrder)
    cells = DataFrame(x=Int64[], y=Int64[], z=Int64[],
    stimStart=Int64[], stimStop=Int64[])
    for (i, g) in enumerate(trialOrder)
        group = target_groups[g]
        start = stimStartIdx[i]
        endT = stimEndIdx[i]
        for (x, y, z) in eachrow(group)
            xR, yR = Int(round(x, digits=0)), Int(round(y, digits=0))
            push!(cells, (xR, yR, z, start, endT))
        end
    end
    cells
end

f_lookup_cellidx(xyzToIdx) = (x,y,z) -> map((a,b,c)->xyzToIdx[(a,b,c)], x, y, z)

"Add a column with index per unique (x,y,z)."
function addCellIdx(cells::DataFrame, xyzToIdx)
    # number each unique (x,y,z) as 1, 2, ... in order observed
    newCells = transform(cells, [:x,:y,:z] => f_lookup_cellidx(xyzToIdx) => :cellIdx)
    # not guranteed to match binary strings used for construction...
    # transform(newCells, :cellIdx => (idx -> num2String(idx)) => :cellStr)
    newCells
end

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

aa2a(array) = copy(hcat(array...)')

function findIdxOfClosestElem(elem, array)
    dist = @. (array - elem)^2
    return sortperm(dist)[1]
end

function getPlaneETLvals(tseries_xml)
    etlVals = sort(unique(round.(parse.(Float64,
        tseries_xml[xpath"//PVStateValue//SubindexedValue[@description='ETL']/@value"]), digits=1)))    
end

"""Convert target_groups z from meters to the plane Int from Imaging.

zOffset is subtracted off. eg holes burn 45um above, so +45*1e-6 zOffset in .mat
and pass 45 to this function.
"""
function mapTargetGroupsToPlane(target_groups, etlVals; is1024=true, zOffset=0.)
    newTargetGroups = deepcopy(target_groups)
    for g in 1:size(target_groups, 1)
        for i in 1:size(target_groups[g], 1)
        # in sean's code, units in meters, so *1e6
            adjustedZ = newTargetGroups[g][i,3] * 1e6 - zOffset
            idx = findIdxOfClosestElem(adjustedZ, etlVals)
            newTargetGroups[g][i,3] = idx
        end
        if is1024
            newTargetGroups[g][:,1:2] *= 2
        else
            @info "will alias fractional target center to nearest Int"
            newTargetGroups[g] .= round.(newTargetGroups[g], digits=0)
        end
    end
    map(g -> Int64.(g), newTargetGroups)
end
    
function getROItrace()
    idx(CI, t) = CartesianIndex(CI.I[1], CI.I[2], t)
    end

function zeroAdjust(image::Array{T,4}) where T <: Real
    # PMT offset
    # @assert sum(image.==0) == 0
    im = image .- 8192
    im[im .< 0] .= 0
    # every other line, reverse column for roundtrip correction
    # for B115, need to start with odd line else fliplr
    im[:,1:2:end,:,:] .= im[:,1:2:end,end:-1:1,:]
    im
end

idxWithTime(CI, t) = CartesianIndex(CI.I[1], CI.I[2], t)

function imageCorrWithMask(images::Array{T,4}, mask; lag=0) where T <: Real
    corImage = zeros(size(images)[1:3]...)
    trace = extractTrace(images, mask)[1:end - lag]
    p = Progress(prod(size(corImage)), 1, "Calc image corr: ")
    @threads for i in CartesianIndices(corImage)
        # corImage[i] = cor(trace, imageJkalmanFilter(view(images,i,:))[31:end])
        corImage[i] = cor(trace, view(images, i, :)[lag + 1:end])
        next!(p)
    end
    corImage
end

function Green(x)
    im = RGB(x)
    RGB(zero(x), im.g, zero(x))
end

    export read_microns_per_pixel,
    read_mask,
    zbrain_units,
    antsApplyTransforms,
    getMaskNameIncludes,
    read_zseries,
    create_nuclei_mask,
    findNeurons,
    create_slm_stim,
    mask2ind,
    readZseriesTiffDir,
    cartesianIndToArray,
    arrCart2arr,
    findNeuronsAaron,
    cartIdx2SeanTarget,
    cartIdx2Array,
    stripLeadingSpace,
    findTTLStarts,
    findTTLEnds,
    readTseriesTiffDir,
    getFramePlane,
    tseriesTiffDirMetadata,
    addTargetsToImage,
    getStimTimesFromVoltages,
    getTimeFromFilename,
    getStimFrameIdx,
    getExpData,
    getTrialOrder,
    loadTseries,
    write_experiment_to_tyh5,
    init_workers,
    trialAverage,
    colorImage,
    countPulsesMaxGap,
    findMatGroups,
    allWithin,
    rectangle,
    extractTrace,
    imageJkalmanFilter,
    stonerStimGroups,
    rpadArray,
    calc_concurrency_score,
    randomSwaps,
    minMaxIntersection,
    count_concurrency,
    noDups,
    printGroups,
    aa2a,
    calcGroupsPerCell,
    stonerPerm,
    makeCellsDF,
    findIdxOfClosestElem,
    mapTargetGroupsToPlane,
    constructGroupMasks,
    imageCorrWithMask,
    Green,
    sym_adjust_gamma,
    zeroAdjust,
    perm, 
    cartIdx2SeanTarget512,
    reconstructA,
    addCellIdx,
    f_lookup_cellidx,
    entangledInfluenceMaps,
    df_f
    # , segment_nuclei
end