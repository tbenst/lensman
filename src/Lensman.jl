module Lensman

using AxisArrays, ANTsRegistration, NIfTI, ImageMagick, Images,
    ImageDraw, ImageFiltering, PyCall, MAT, Dates, DataStructures,
    Statistics, SharedArrays, CSV, DataFrames
import Base.Threads.@threads
using Distributed
# https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules

peak_local_max = PyNULL()
disk = PyNULL()
match_template = PyNULL()
np = PyNULL()
autolevel_percentile = PyNULL()
scipy = PyNULL()

function __init__()
    pyimport("skimage")
    copy!(peak_local_max, pyimport("skimage.feature").peak_local_max)
    copy!(disk, pyimport("skimage.morphology").disk)
    copy!(match_template, pyimport("skimage.feature").match_template)
    copy!(np, pyimport("numpy"))
    copy!(scipy, pyimport("scipy"))
    pyimport("scipy.ndimage")
    copy!(autolevel_percentile, pyimport("skimage.filters.rank").autolevel_percentile)
end
    

include("Bruker.jl")
# include("segment.jl")

zbrain_units = (0.798μm, 0.798μm, 2μm)
tmppath = ANTsRegistration.userpath()

"Read (sparse) mask from Zbrain atlas and reshape to dense Array."
function read_mask(masks, idx; W=621, H=1406, Z=138, units=zbrain_units)
    example_mask = reshape(Array(masks["MaskDatabase"][:,idx]), (H, W, Z))
    # reshape as [621, 1406, 138]
    example_mask = permutedims(example_mask, (2,1,3))
    AxisArray(example_mask, (:y, :x, :z), units)
end

"Combine Zbrain masks whose name includes the given string."
function getMaskNameIncludes(masks, maskNameIncludes, units=zbrain_units)
    masks_raphe_idxs = findall(occursin.(maskNameIncludes, masks["MaskDatabaseNames"])) |>
        x -> map(ci -> ci[2], x)
    mask = foldl((a,n)->a .| read_mask(masks, n), masks_raphe_idxs;
        init=read_mask(masks,masks_raphe_idxs[1]))
    mask = AxisArray(Float32.(mask), (:y, :x, :z), zbrain_units)
end

function antsApplyTransforms(fixedPath::String, moving::AxisArray, transformPath::String)
    movingPath = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
    antsApplyTransforms(fixedPath, movingPath, transformPath)
end

function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath::String)
    maskoutname = joinpath(tmppath, ANTsRegistration.randstring(10)*".nii.gz")
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
    @threads for z in 1:size(tiff_files,1)
        zseries[:,:,z] = ImageMagick.load(tiff_files[z])
    end
    return zseries
end

function create_nuclei_mask(microscope_units)
    nucleus_r = (2.75μm/microscope_units[1])μm
    space = 0.3μm # WARN: need to adjust CirclePointRadius, too
    cell_r = nucleus_r + space
    template_size = Int32(ceil(2*(cell_r)/microscope_units[1]))
    # ensure odd size
    template_size += abs((template_size) % 2 - 1)
    # we make large and later resize for anti-aliasing
    nucleus_template = Gray.(zeros(513,513))
    draw!(nucleus_template, Ellipse(CirclePointRadius(256,256,192)))
    nucleus_template = Float64.(nucleus_template)
    # nucleus_template[nucleus_template.==0.] .= -1.
    sz = (template_size, template_size)
    σ = map((o,n)->0.75*o/n, size(nucleus_template), sz)
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

    stim_struct = OrderedDict{String, Any}()
    # field order: im, targets, zrange, targetsGroup
    exp_key = "experimentS"
    # exp_key = "etb"
    stim_struct[exp_key] = Dict{String, Any}()
    out_mat = stim_struct[exp_key]

    out_mat["im"] = Array{Any}(nothing,1,2)
    # H, W = 1024, 1024
    H, W = 512, 512
    out_mat["im"][1] = zeros(H,W)
    out_mat["im"][2] = zeros(H,W)
    @assert all(targets .<512)

    out_mat["targets"] = Array{Any}(nothing,1,1)
    out_mat["targets"][1] = copy(targets) # ensure not an Adjoint type
    out_mat["zrange"] = [minimum(targets[:,3]) maximum(targets[:,3])]
    # channel 1 (red) and 2 (green), used for GUI preview

    ntargets = size(targets,1)
    
    # perhaps this is optional...?
    # targetsGroup alternates between SLM1&2 (eg Float32(N,3) and zeros(0,0)) and each one is 2ms apart (4ms between first two SLM1 patterns)
    # numGroups = 99 # not sure why this varies...
    numGroups = 4
    out_mat["targetsGroup"] = Array{Any}(nothing,1,numGroups)
    out_mat["targetsGroup"][1] = copy(targets) # stim at +0ms

    for i in 2:2:numGroups
        # SLM2 is not used
        out_mat["targetsGroup"][i] = zeros(0,0)
    end
    for i in 3:2:numGroups
        # future stim progressions go here
        out_mat["targetsGroup"][i] = zeros(0,0)
    end

    # println("construct 10Hz pattern ")
    # out_mat["targetsGroup"][50] = copy(targets) # stim at +100ms

    matwrite(matpath, stim_struct)
    # strip /mnt/deissero... this is needed since Matlab on windows will read the file
    return "O:\\"*replace(matpath[14:end], "/" => "\\")
end

# ::Array{String,1}
function create_trials_txt(targets_mats, outname::String;
        outdir = "/mnt/deissero/users/tyler/slm/masks/")
    today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
    today_dir = mkpath(joinpath(outdir, today_str))
    txtpath = joinpath(today_dir, "$outname.txt")
    
    txt_io = open(txtpath, "w")
    for matpath in targets_mats
        println(txt_io, matpath*"\t1")
    end
    close(txt_io)
    return txtpath
end
"Outname is for local filesystem (build platform)"
function create_slm_stim(target_groups, outname::String)
    targets_mats = []
    for (i, targets) in enumerate(target_groups)
        name = "$(outname)_group_$i"
        push!(targets_mats,create_targets_mat(targets, name))
    end
    trials_txt = create_trials_txt(targets_mats, outname)
    return targets_mats, trials_txt
end


function cartIdx2Array(x::Vector{<:CartesianIndex})
    # make tall matrix
    permutedims(reduce(hcat, collect.(Tuple.(x))), (2,1))
end

cartesianIndToArray = cartIdx2Array

function readZseriesTiffDir(tifdir; contains="Ch3")
    tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") & occursin(contains, x),
        readdir(tifdir)))
    tif0 = ImageMagick.load(tiff_files[1])
    H, W = size(tif0)

    zseries = zeros(Normed{UInt16,16}, H, W, size(tiff_files, 1))
    @threads for z in 1:size(tiff_files,1)
        zseries[:,:,z] = ImageMagick.load(tiff_files[z])
    end
    zseries
end

"From bruker tif path, return (frameNum, planeNum)"
function getFramePlane(tifPath)
    m = match(r".*Cycle(\d+)_Ch\d_(\d+).ome.tif", tifPath)
    frameNum, planeNum = parse(Int,m[1]), parse(Int,m[2])
end


"""Prepare to read data in parallel.

Usage:
H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(".")
addprocs(18)
@everywhere using Lensman
tseries = SharedArray{Normed{UInt16,16},4}((H, W, Z, T),
    init = S -> _readTseriesTiffDir(S, framePlane2tiffPath, Z, T))

    """
function tseriesTiffDirMetadata(tifdir, containsStr = "Ch3")
    tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") &
        occursin(containsStr, x), readdir(tifdir)))
    framePlane2tiffPath = Dict(getFramePlane(tp) => tp for tp in tiff_files)
    framePlanes = hcat(collect.(keys(framePlane2tiffPath))...)
    T = maximum(framePlanes[1,:])
    Z = maximum(framePlanes[2,:])
    if T==1
        # single plane imaging, so T is last int in filename...
        framePlane2tiffPath = Dict((k[2],k[1]) => v for (k,v) in framePlane2tiffPath)
        T, Z = Z, T
    end
    tif0 = ImageMagick.load(tiff_files[1])
    H, W = size(tif0)
    H, W, Z, T, framePlane2tiffPath
end

function _readTseriesTiffDir(S::SharedArray, framePlane2tiffPath, Z, T)
    nProcs = size(procs(S),1)
    i = indexpids(S)
    myT = collect(i:nProcs:T)
    println("process $i will open $(size(myT,1)) tiffs")
    for t in myT
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            S[:,:,z,t] .= ImageMagick.load(tp)
        end
    end
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
    frameStarted = diff(voltages .> std(voltages).*3)
    frameStartIdx = findall(frameStarted .== 1) .+ 1
end

function findTTLEnds(voltages)
    frameStarted = diff(voltages .< std(voltages).*3)
    frameStartIdx = findall(frameStarted .== 1) .+ 1
end

"""Pass image (HxW) and targetLocs (Nx2), and return image HxW"

targetSize=4.41 comes from 7um spiral on 25x objective @ 1x (0.63 μm/px)
"""
function addTargetsToImage(image::Array{<:Real,2}, targetLocs::Array{<:Real,2};
    targetSize::Real=11.11, alpha=0.4)
    im = RGB.(image)
    stim = zeros(Gray,size(image))
    for (x,y) in eachrow(targetLocs)
        draw!(stim, Ellipse(CirclePointRadius(y,x,targetSize)))
    end
    channelview(im)[[1,3],:,:,:] .= 0
    channelview(im)[1,:,:,:] .= stim*alpha
    im
end


"Convert mask to an Array of indices, similar to argwhere in Python."
function mask2ind(mask)
    transpose(hcat(collect.(Tuple.(findall(mask)))...))
end

# /2: account for sean's code using 512 x 512 space
# while we image in 1024x1024
cartIdx2SeanTarget(ind, z) = hcat(reverse(collect(Tuple(ind))/2)..., z)

""
function getStimTimesFromVoltages(voltageFile, Z::Int)
    voltages = CSV.File(open(read, voltageFile)) |> DataFrame
    rename!(voltages, [c => stripLeadingSpace(c) for c in names(voltages)])

    # plot(voltages["Time(ms)"][1:10000], voltages["frame starts"][1:10000])
    frameStarted = diff(voltages[!,"frame starts"] .> std(voltages[!,"frame starts"]).*3)
    frameStartIdx = findall(frameStarted .== 1) .+ 1

    @warn "hardcodes 5vol stim (1s @ 5Hz)"
    stimDur = 5
    stimStartIdx = findTTLStarts(voltages[!,"respir"])[stimDur:stimDur:end]
    stimEndIdx = findTTLEnds(voltages[!,"respir"])[stimDur:stimDur:end]
    stimStartFrameIdx = [Int.(floor.(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimStartIdx]
    stimEndFrameIdx = [Int.(ceil(searchsortedfirst(frameStartIdx, s)/Z)) for s in stimEndIdx]
    stimStartFrameIdx, stimEndFrameIdx
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
    getStimTimesFromVoltages
    # , segment_nuclei
end