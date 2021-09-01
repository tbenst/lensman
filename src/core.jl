peak_local_max = PyNULL()
disk = PyNULL()
match_template = PyNULL()
np = PyNULL()
autolevel_percentile = PyNULL()
scipy = PyNULL()
machinery = PyNULL()
py_utils = PyNULL()


function __init__()
    pyimport_conda("skimage", "scikit-image", "tables")
    # pyimport("skimage")
    copy!(peak_local_max, pyimport("skimage.feature").peak_local_max)
    copy!(disk, pyimport("skimage.morphology").disk)
    copy!(match_template, pyimport("skimage.feature").match_template)
    copy!(np, pyimport("numpy"))
    copy!(scipy, pyimport("scipy"))
    pyimport("scipy.ndimage")
    copy!(autolevel_percentile, pyimport("skimage.filters.rank").autolevel_percentile)

    # Nifty util to load local Python source code
    copy!(machinery, pyimport("importlib.machinery"))
    loader = machinery.SourceFileLoader("py_utils", joinpath(dirname(@__FILE__), "py_utils.py"))
    copy!(py_utils, loader.load_module("py_utils"))
end


zbrain_units = (0.798μm, 0.798μm, 2μm)

try
    global tmppath = ANTsRegistration.userpath()
catch
    @warn "ANTs does not officially support Windows. Disabling."
    global tmppath = mktemp()
end

function zbrain_vec2mat(sparse_vec; W=621, H=1406, Z=138, units=zbrain_units,
        rostral=:left, dorsal=:down)
    mask = reshape(Array(sparse_vec), (H, W, Z))
    # reshape as [621, 1406, 138]
    mask = permutedims(mask, (2, 1, 3))
    # zbrain is facing left and bottom to top, meaning idx 1 is deepest plane
    if rostral == :right
        mask = reverse(mask,dims=2);
    end
    if dorsal == :up
        mask = reverse(mask,dims=3);
    end
    AxisArray(mask, (:y, :x, :z), units)
end

"Read (sparse) mask from Zbrain atlas and reshape to dense Array."
function read_mask(masks, idx::Int; W=621, H=1406, Z=138, units=zbrain_units,
        rostral=:left, dorsal=:down)
    zbrain_vec2mat(masks["MaskDatabase"][:,idx];
        W=W, H=H, Z=Z, units=units, rostral=rostral,dorsal=dorsal)
end

function read_mask(masks, names::AbstractString; kwargs...)
    mask_names = masks["MaskDatabaseNames"]
    # we change name when resaving h5 file
    new_mask_names = replace.(mask_names, "/" => "_")[1,:]
    # lets resolve for either
    is_name = mask_names .== names
    new_is_name = new_mask_names .== names
    if sum(is_name) == 1
        idx = findall(is_name)[1]
    elseif sum(new_is_name) == 1
        idx = findall(new_is_name)[1]
    else
        error("could not find (unique) name, $name")
    end

    read_mask(masks, idx; kwargs...)
end

"Combine Zbrain masks whose name includes the given string."
function getMaskNameIncludes(masks, maskNameIncludes, units=zbrain_units)
    masks_raphe_idxs = findall(occursin.(maskNameIncludes, masks["MaskDatabaseNames"])) |>
        x -> map(ci -> ci[2], x)
    mask = foldl((a, n) -> a .| read_mask(masks, n), masks_raphe_idxs;
        init=read_mask(masks, masks_raphe_idxs[1]))
    mask = AxisArray(Float32.(mask), (:y, :x, :z), zbrain_units)
end

"""Apply one transform, and write moving to a tmp file.

5 different methods, depending on if 1 transform (rigid only) or 2
transforms (rigid + nonrigid), and if passing a filepath or an array.
"""
function antsApplyTransforms(fixedPath::String, moving::AxisArray, transformPath::String)
    movingPath = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
    result = antsApplyTransforms(fixedPath, movingPath, transformPath)
    rm(movingPath)
    result
end

"Apply one transform."
function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath::String)
    maskoutname = joinpath(tmppath, ANTsRegistration.randstring(10) * ".nii.gz")
    result = antsApplyTransforms(fixedPath, movingPath, transformPath, maskoutname)
    rm(maskoutname)
    result
end

# function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath1::String, transformPath2::String)
#     maskoutname = joinpath(tmppath, ANTsRegistration.randstring(10) * ".nii.gz")
#     antsApplyTransforms(fixedPath, movingPath, transformPath1, transformPath2, maskoutname)
# end

"Apply one transform and save file with given outname."
function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath::String, maskoutname::String,
        ants_path = "/opt/ANTs/install/bin/antsApplyTransforms")
    run(`$ants_path --float -d 3 -i $movingPath -r $fixedPath -t $transformPath -o $maskoutname`)
    niread(maskoutname)
end

"Apply two transforms, writing AxisArrays to tmp files."
function antsApplyTransforms(fixed::AxisArray, moving::AxisArray,
        transformPath1::String, transformPath2::String,
        maskoutname = joinpath(tmppath, ANTsRegistration.randstring(10) * ".nii.gz"),
        ants_path = "/opt/ANTs/install/bin/antsApplyTransforms")
    fixedPath = joinpath(tmppath, ANTsRegistration.write_nrrd(fixed))
    movingPath = joinpath(tmppath, ANTsRegistration.write_nrrd(moving))
    result = antsApplyTransforms(fixedPath, movingPath, transformPath1, transformPath2,
        maskoutname, ants_path)
    rm(fixedPath)
    rm(movingPath)
    result
end

"Apply two transforms on existing fixed & moving files."
function antsApplyTransforms(fixedPath::String, movingPath::String, transformPath1::String, transformPath2::String, maskoutname::String,
        ants_path = "/opt/ANTs/install/bin/antsApplyTransforms")
    run(`$ants_path --float -d 3 -i $movingPath -r $fixedPath -t $transformPath2 -t $transformPath1 -o $maskoutname`)
    niread(maskoutname)
end

"Tranform all 4 channels from e.g. Olympus."
function multimap_transforms(fixed, oir_img, mm_affine_transform, mm_SyN_transform)
    registered = zeros(size(fixed)[1:2]...,4,size(fixed,3))
    # TODO: perhaps could cache this...
    @threads for c in 1:4
        registered[:,:,c,:] = antsApplyTransforms(fixed, oir_img[:,:,c,1,:],
            mm_affine_transform, mm_SyN_transform)
    end
    registered
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
function findNeurons(plane; thresh_adjust=1.2, featSize=2, maxiSize=4)
    # mask out background
    adj_plane = adjust_histogram(imadjustintensity(plane), GammaCorrection(0.5))
    threshold = otsu_threshold(adj_plane) * thresh_adjust
    otsu_mask = adj_plane .> threshold
    img = adj_plane .* otsu_mask
    feat_img = match_template(Float64.(img), disk(featSize), pad_input=true, mode="constant", constant_values=0)
    local_maxi = peak_local_max(feat_img, indices=false, footprint=disk(maxiSize), exclude_border=false, threshold_rel=.1)
end

function cartIdx2Array(x::Vector{<:CartesianIndex})
    # make tall matrix
    permutedims(reduce(hcat, collect.(Tuple.(x))), (2, 1))
end

"ImageJ 'kalman filter'"
function imageJkalmanFilter(ys::Array{T,1};gain=0.8,useBurnIn=true) where T
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

"TODO: improve memory efficiency."
function imageJkalmanFilter(ndarray;gain=0.8,useBurnIn=true) where T
    denoised = zeros(size(ndarray)...)
    @threads for i in CartesianIndices(denoised[:,:,:,1])
        denoised[i,:] = imageJkalmanFilter(ndarray[i,:])
    end
    denoised
end
    
function extractTrace(tseries, mask::Array{Bool,T}) where T
    extractTrace(tseries, findall(mask))
end

function extractTrace(tseries, cartIdxs::Array{T,1}) where T <: CartesianIndex
    nT = size(tseries, ndims(tseries))
    trace = zeros(Float64, nT)
    for t in 1:nT
        # trace[t] = mean(reinterpret(UInt16, view(tseries, cartIdxs, t)))
        # trace[t] = mean(view(tseries, cartIdxs, t))
        trace[t] = mean(selectdim(tseries, 4, t)[cartIdxs])
    end
    trace
end

# "Calculate df/f, assuming last dimension is time."
# function df_f(f_timeseries, f0_timeseries)
#     d = ndims(f_timeseries)
#     f = selectdim(mean(f_timeseries, dims=d), d, 1)
#     f0 = selectdim(mean(f0_timeseries, dims=d), d, 1)
#     @. (f - f0) / f0
# end

"Calculate df/f, assuming last dimension is time."
function timeseries_df_f(f_timeseries)
    d = ndims(f_timeseries)
    f = f_timeseries
    f0 = minimum(f_timeseries, dims=d)
    @. (f - f0) / f0
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


stripLeadingSpace(s) = (s[1] == ' ') ? s[2:end] : s

function findTTLStarts(voltages)
    frameStarted = diff(voltages .> std(voltages) .* 3)
    # add one due to diff
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

function getExpData(tseries_path::String)
    open(xmlPath, "r") do io
        tseries_xml = read(io, String)
        xp_parse(tseries_xml)
    end
end
function getExpData(tseries_xml::LibExpat.ETree)
    # warning: find has different syntax from xpath!!!
    expDate = LibExpat.find(tseries_xml, "/PVScan[1]{date}")
    expDate = DateTime(expDate, "m/d/y I:M:S p")
    # much faster by prepending /PVScan/PVStateShard[1] to search first frame only
    framePeriod = parse(Float64, tseries_xml[xpath"/PVScan/PVStateShard[1]/PVStateValue[@key='framePeriod']/@value"][1])
    frameRate = 1 / framePeriod
    etlVals = nothing
    try
        etlVals = getPlaneETLvals(tseries_xml)
    catch
        # TODO: should print exception
        etlVals = nothing
    end
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
function getTrialOrder(slm_dir, expDate)
    date_str = Dates.format(expDate, "dd-u-Y")
    slmExpDir = joinpath(slm_dir,date_str)
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
    if (timeDiff > Second(120)) | (timeDiff < Second(0))
        @warn "started SLM $(timeDiff) seconds before imaging on $date_str. Is this correct?"
    end
    trialOrderTxt = trialOrderTxts[dtIdx]

    trialOrder = CSV.File(trialOrderTxt, delim="\t", header=true) |> Tables.matrix
    @assert size(trialOrder, 1) == 1
    trialOrder[1,:], joinpath(splitpath(trialOrderTxt)[1:end - 1]...)
end

# function _loadTSeries!(tseries, t, framePlane2tiffPath, p, Z)
#     for z in 1:Z
#         tp = framePlane2tiffPath[t,z]
#         tseries[:,:,z,t] .= reinterpret(UInt16, ImageMagick.load(tp))
#     end
#     next!(p)
# end

# "Read a 4D (HWZT) tiff stack in parallel from a folder."
# function loadTseries(tifdir, containsStr::String="Ch3")
#     H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
#     tseries = Array{UInt16}(undef, H, W, Z, T)
#     memory_size_bytes = prod(size(tseries)) * 2
#     memory_size_gb = round(memory_size_bytes / 1024^3, digits=1)
#     println("estimated memory usage: $memory_size_gb")
#     p = Progress(T, 1, "Load Tseries: ")
#     Threads.@sync for t in 1:T
#         Threads.@spawn _loadTSeries!(tseries, t, framePlane2tiffPath, p, Z)
#     end
# end

"Read a 4D (HWZT) tiff stack in parallel from a folder."
function loadTseries(tifdir, containsStr::String="Ch3"; missing=false)
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
    if missing
        tseries = Array{Union{Missing,UInt16}}(missing, H, W, Z, T)
    else
        tseries = Array{UInt16}(undef, H, W, Z, T)
    end
    memory_size_bytes = prod(size(tseries)) * 2
    memory_size_gb = round(memory_size_bytes / 1024^3, digits=1)
    println("estimated memory usage: $memory_size_gb")
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

"Sparse, memory efficent average of Tseries"
function avgTseries(tifdir, containsStr::String="Ch3"; every=100)
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
    newT = Int(ceil(T / every))
    tseries = Array{UInt16}(undef, H, W, Z, newT)
    memory_size_bytes = prod(size(tseries)) * 2
    memory_size_gb = round(memory_size_bytes / 1024^3, digits=1)
    println("estimated memory usage: $memory_size_gb")
    p = Progress(T, 1, "Average tseries:")
    times = collect(1:every:T)
    # @threads for t in 1:every:T
    @threads for i in 1:length(times)
        t = times[i]
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            tseries[:,:,z,i] .= reinterpret(UInt16, ImageMagick.load(tp))
        end
        next!(p)
    end
    (sum(tseries, dims=4) ./ newT)[:,:,:,1]
end

"Load brightness over time of mask"
function loadBOT(tifdir; freduce=mean, containsStr::String="Ch3")
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tifdir, containsStr)
    ex_im = ImageMagick.load(framePlane2tiffPath[1,1])
    ex_reduced = freduce(ex_im)
    tseries = Array{eltype(ex_reduced)}(undef, size(ex_reduced)..., Z, T)
    p = Progress(T, 1, "Load Tseries: ")

    @threads for t in 1:T
        for z in 1:Z
            tp = framePlane2tiffPath[t,z]
            tseries[..,z,t] .= freduce(ImageMagick.load(tp))
        end
        next!(p)
    end
    tseries
end

function get_slm_stim_masks(tif_dir, slm_dir)
    """
    # Arguments
    - `tif_dir::String`: path to directory holding experiment data (<fish_dir>/<exp_name>)
    - `slm_dir::String`: path to directory containing SLM metadata for experiments
       (e.g. "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment")

    # note
    - This function is currently only valid for use with the 25x objective and NOT with
      the 16x (14.4x) objective
    """
    @error "this appears to give the wrong results... -Tyler"
    return (nothing, nothing)
    # we avoid using tseries since so slow to read in, and don't want to error on
    # directory not mounted, etc.
    H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tif_dir)

    # TODO(allan.raventos): really should get rid of this `is_1024`
    if W == 512
        is_1024 = false
    elseif W == 1024
        is_1024 = true
    else
        error("Only supporting (H, W) in {(512, 512), (1024, 1024)}")
    end

    @warn "lateral unit hard coded for B115 1x zoom"
    lateral_unit = microscope_lateral_unit(W)

    # Parse stimulation start and end times from voltage recordings
    voltage_file = glob("*VoltageRecording*.csv", tif_dir)[1]
    stim_start_idx, stim_end_idx = getStimTimesFromVoltages(voltage_file, Z)
    allWithin(diff(stim_start_idx),0.05)

    # If empty we have no stimulations so return nothing
    if length(stim_start_idx) == 0
        return (nothing, nothing)
    end

    # Read SLM stim files
    data_folders = splitpath(tif_dir)
    # xml_path is e.g. <fish_dir>/<exp_name>/<exp_name>.xml
    xml_path = joinpath(data_folders..., data_folders[end] * ".xml")
    exp_date, frame_rate, etl_vals = getExpData(xml_path)

    # If imaging many planes, may need to read another xml file since didn't get all planes in first file
    # Volume gets imaged at 1/Zth of the frame rate
    @assert length(etl_vals) == Z
    @info "assume imaging from top-down"

    # Load SLM files (TODO: handle this constant path better)
    slm_exp_dir = joinpath(slm_dir, Dates.format(exp_date, "dd-u-Y"))
    trial_order, slm_exp_dir = getTrialOrder(slm_exp_dir, exp_date)
    n_stimuli = maximum(trial_order)
    n_trials = size(trial_order, 1)
    @assert n_trials == size(stim_start_idx, 1)
    target_radius_px = spiral_size(exp_date, lateral_unit)


    mat = matread.(findMatGroups(slm_exp_dir)[1])
    max_targets_center = maximum(map(mat->maximum(maximum.(matread(mat)["cfg"]["exp"]["maskS"]["targetsCenter"])), findMatGroups(slm_exp_dir)))
    @assert  max_targets_center < 512 "bad targets center value, don't trust results..."
    slm_num = getSLMnum(mat)
    z_offset = getzoffset(exp_date, slm_num)

    firstTargetGroup = matread.(findMatGroups(slm_exp_dir)[1])
    powerPerCell = firstTargetGroup["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"]
    slm1Power, slm2Power = slmpower(exp_date)
    if slm_num == 1
        slmpowerPerCell = slm1Power * powerPerCell / 1000
    elseif slm_num == 2
        slmpowerPerCell = slm2Power * powerPerCell / 1000
    end

    

    # Load target groups
    # Assumes no sequence of stim
    # for 5Hz clock

    target_groups = []
    group_stim_freq = []
    for mat in matread.(findMatGroups(slm_exp_dir))
        push!(target_groups, mat["cfg"]["maskS"]["targets"][1])
        push!(group_stim_freq, getMatStimFreq(mat))
    end
        
    @assert length(target_groups) == n_stimuli "trialOrder has max of $n_stimuli, but target_groups is length $(length(target_groups))"
    ## Array of matrices indicating x, y z (or was it y, x, z..?)
    targets_with_plane_index = mapTargetGroupsToPlane(target_groups, etl_vals,
        is1024=is_1024, zOffset=z_offset)
    targets_with_plane_index = map(x->Int.(round.(x, digits=0)), targets_with_plane_index)
    # Generate mask of targeted locations
    stim_masks = constructGroupMasks(targets_with_plane_index, H, W, Z)

    # Create array holding at position i a 0 if no stim was begin applied at timestep i and
    # stim_mask_idx otherwise (meaning that `stim_masks[stim_mask_idx]` was used)
    stim_used_at_each_timestep = zeros(Int64, T)
    for (_start, _end, trial_stim_idx) in zip(stim_start_idx, stim_end_idx, trial_order)
        for idx in _start:_end
            stim_used_at_each_timestep[idx] = trial_stim_idx
        end
    end
    stim_masks, stim_used_at_each_timestep;
end


function write_experiment_to_tyh5(tseries, stim_masks, stim_used_at_each_timestep, output_path; compression_level=3)
    """
    # Arguments
    - `tseries::Array{UInt16, 4}`: (H, W, Z, T) fluorescence volume time series
    - `stim_masks`
    - `stim_used_at_each_timestep`
    - `output_path::String`: full path to write .ty.h5 file to
    - `compression_level::Int`: compression level to use for blosc:zstd when writing out data
    """
    py_utils.write_experiment_to_tyh5(
        tseries,
        output_path,
        stim_masks=stim_masks,
        stim_used_at_each_timestep=stim_used_at_each_timestep,
        compression_level=compression_level
    )

end

"Calculate average response for each unique stimuli"
function trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder; pre=16, post=16)
    nStimuli = maximum(trialOrder)
    nTrials = size(trialOrder, 1)
    nTrialsPerStimulus = [sum(trialOrder .== i) for i in 1:nStimuli]
    @assert nTrials == size(stimStartIdx, 1) # check TTL start-times match

    ## avg stim effec
    trialTime = minimum(stimEndIdx .- stimStartIdx) + pre + post + 1
    # HWZCT
    avgStim = zeros(size(tseries)[1:3]..., nStimuli, trialTime);
    p = Progress(length(stimStartIdx), 1, "trial average:")
    # threads crashes HDF5
    # @threads for i in 1:length(stimStartIdx)
    for i in 1:length(stimStartIdx)
        start = stimStartIdx[i]
        stop = start - pre + trialTime - 1
        trialType = trialOrder[i]
        avgStim[:,:,:,trialType,:] .+= tseries[:,:,:,start - pre:stop]
        next!(p)
    end
    for i in 1:nStimuli
        avgStim[:,:,:,i,:] ./= nTrialsPerStimulus[i]
    end
    avgStim
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
function makeCellsDF(target_groups, stimStartIdx, stimEndIdx, trialOrder, group_stim_freq)
    cells = DataFrame(x=Int64[], y=Int64[], z=Int64[], stimGroup=Int64[],
        stimStart=Int64[], stimStop=Int64[], cellID=Int64[], stimNum=Int64[])
    cellIDs = Dict()
    for (i, g) in enumerate(trialOrder)
        group = target_groups[g]
        start = stimStartIdx[i]
        endT = stimEndIdx[i]
        for (x, y, z) in eachrow(group)
            if (x,y,z) in keys(cellIDs)
                cellID = cellIDs[(x,y,z)]
            else
                cellID = length(cellIDs) + 1
                cellIDs[(x,y,z)] = cellID
            end
            xR, yR = Int(round(x, digits=0)), Int(round(y, digits=0))
            push!(cells, (xR, yR, z, i, start, endT, cellID, g))
        end
    end
    cells[!, :stimFreq] = map(g->group_stim_freq[trialOrder][g], cells.stimGroup)
    cells
end


"""Create DataFrame with row per individual cell stimulated.

Used for recordings with many trials in Prarie View.
"""
function makeCellsDF(target_groups, stimStartIdx::Int, stimEndIdx::Int, trialOrder)
    cells = DataFrame(x=Int64[], y=Int64[], z=Int64[], trialNum=Int64[],
        stimStart=Int64[], stimStop=Int64[], cellID=Int64[])
    cellIDs = Dict()
    for (i, g) in enumerate(trialOrder)
        group = target_groups[g]
        for (x, y, z) in eachrow(group)
            if (x,y,z) in keys(cellIDs)
                cellID = cellIDs[(x,y,z)]
            else
                cellID = length(cellIDs) + 1
                cellIDs[(x,y,z)] = cellID
            end
            xR, yR = Int(round(x, digits=0)), Int(round(y, digits=0))
            push!(cells, (xR, yR, z, i, stimStartIdx, stimEndIdx, cellID))
        end
    end
    cells
end


"""Add stim-evoked df/f to cells DataFrame.

padding: ignore +/- that amount around time of stim
"""
function add_df_f_to_cells(tseries, cells, roiMask;
    winSize=3, padding=2, threads=false, ϵ=1.0
)
    cellDf_f = Float64[]
    cell_f = Float64[]
    cell_f0 = Float64[]
    cell_area = Float64[]
    stim_grouped_df = groupby(cells, :stimStart)
    # we use to calculate the row in cellDf_f
    lengths = combine(stim_grouped_df, :x => length).x_length
    T = size(tseries, ndims(tseries))
    # kalman_burnin = 5
    # we don't use @threads in case tseries isn't thread-safe
    p = Progress(length(stim_grouped_df); desc="df/f:")
    # @showprogress for (g, grouped_df) in enumerate(stim_grouped_df)
    if threads
        mapper = Folds.map
    else
        mapper = map
    end
    mapper(1:length(stim_grouped_df)) do g
        grouped_df = stim_grouped_df[g]
        stimStart, stimStop = grouped_df[1, [:stimStart, :stimStop,]]
        # theStart = stimStart - winSize - kalman_burnin - padding
        theStart = stimStart - winSize - padding
        @assert theStart >= 1
        theEnd = stimStop + winSize + padding
        @assert theEnd <= T
        duration = theEnd - theStart
        # indexing can be slow for some types, so do in outer loop
        tseries_trial = tseries[:, :, :, theStart:theEnd]
        for i in 1:size(grouped_df, 1)
            x, y, z, cellID = grouped_df[i, [:x, :y, :z, :cellID]]
            mask = roiMask[cellID]
            neuronTrace = extractTrace(tseries_trial, mask)
            # neuronTrace = imageJkalmanFilter(neuronTrace)[kalman_burnin+1:end]
            f = mean(neuronTrace[end-winSize+1:end])
            
            # problem: if neuron happens to be active and declining before stim,
            # we may miss
            f0 = mean(neuronTrace[1:winSize])

            # solution: use percentile f0
            # f0 = percentile(neuronTrace,0.1)
            df_f = (f - f0) / (f0 + ϵ)
            if isnan(df_f)
                @show f, f0, i, x, y, z, stimStart, stimStop, theStart, stimStart - 1
            end
            area = length(mask)
            push!(cellDf_f, df_f)
            push!(cell_f, f)
            push!(cell_f0, f0)
            push!(cell_area, area)
        end
        next!(p)
        nothing
    end

    # ensure same row order as cellDf_f
    newCells = transform(stim_grouped_df)
    if "df_f" in names(cells)
        newCells.df_f = cellDf_f
        return newCells
    else
        return insertcols!(newCells, size(cells, 2) + 1,
            :df_f => cellDf_f, :f => cell_f, :f0 => cell_f0, :area => cell_area)
    end
end

"""Create binary mask of HxWxZ for each neuron location.
"""
function constructROImasks(cells, H, W, Z, targetSizePx)
    @warn "now must index by cellID"
    nCells = maximum(cells.cellID)
    roiMask = Dict()
    @threads for cell in eachrow(cells)
        (x, y, z, group, stimStart, stimStop, cellID) = cell[[:x,:y,:z,:stimGroup,:stimStart,:stimStop,:cellID]]
        mask = Gray.(zeros(Bool, H, W, Z))
        draw!(view(mask, :, :, z), Ellipse(CirclePointRadius(x, y, targetSizePx)))
        # TODO: should we replace this with sparse array..? Prob should benchmark.
        mask = findall(channelview(mask))
        roiMask[cellID] = mask
    end
    roiMask
end

function constructROImasks(cells::DataFrame, H, W, targetSizePx)
    @warn "now must index by cellID"
    nCells = maximum(cells.cellID)
    roiMask = Dict()
    @threads for cell in eachrow(cells)
        (x,y,z,trialNum,stimStart,stimStop, cellID) = cell[[:x,:y,:z,:trialNum,:stimStart,:stimStop,:cellID]]
        mask = Gray.(zeros(Bool, H,W))
        draw!(view(mask,:,:), Ellipse(CirclePointRadius(x,y,targetSizePx)))
        mask = findall(channelview(mask))
        roiMask[cellID] = mask
    end
    roiMask
end

function constructROImasks(neuron_locs::Vector{CartesianIndex{2}}, H, W, targetSizePx)
    roiMask = []
    for (i,loc) in enumerate(neuron_locs)
        y,x = Tuple(loc)
        mask = Gray.(zeros(Bool, H,W))
        draw!(view(mask,:,:), Ellipse(CirclePointRadius(x,y,targetSizePx)))
        mask = findall(channelview(mask))
        push!(roiMask, mask)
    end
    roiMask
end



f_lookup_cellidx(xyzToIdx) = (x, y, z) -> map((a, b, c) -> xyzToIdx[(a, b, c)], x, y, z)

"Add a column with index per unique (x,y,z)."
function addCellIdx(cells::DataFrame, xyzToIdx)
    # number each unique (x,y,z) as 1, 2, ... in order observed
    newCells = transform(cells, [:x,:y,:z] => f_lookup_cellidx(xyzToIdx) => :cellIdx)
    # not guranteed to match binary strings used for construction...
    # transform(newCells, :cellIdx => (idx -> num2String(idx)) => :cellStr)
    newCells
end

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

"Array of arrays to array."
aa2a(array) = copy(hcat(array...)')
# aa2a(array) = collect(Iterators.flatten(array))

function findIdxOfClosestElem(elem, array)
    dist = @. (array - elem)^2
    return sortperm(dist)[1]
    end

function getPlaneETLvals(tseries_xml)
    etlVals = sort(unique(round.(parse.(Float64,
        tseries_xml[xpath"//PVStateValue//SubindexedValue[@description='ETL']/@value"]), digits=1)))
    if (sum(etlVals .=== 0.0)==1) & (sum(etlVals .=== -0.0)==1)
        return etlVals[(~).(etlVals .=== -0.0)]
    else
        etlVals
    end
end

"""Convert target_groups z from meters to the plane Int from Imaging.

zOffset is subtracted off. eg holes burn 45um above, so +45*1e-6 zOffset in .mat
and pass 45 to this function.
"""
function mapTargetGroupsToPlane(target_groups, etlVals; is1024=true, zOffset=0.)
    # @show typeof(target_groups), typeof(etlVals), typeof(is1024), typeof(z_offset)
    
    # TODO: bizarre bug in Thunks...?
    # totally unknown why this reify call is needed
    # this is a total hack
    @show typeof(etlVals)
    z_offset = reify(zOffset)
    # maybe deepcopy is breaking...?
    newTargetGroups = deepcopy(target_groups)
    for g in 1:size(target_groups, 1)
        for i in 1:size(target_groups[g], 1)
        # in sean's code, units in meters, so *1e6
            adjustedZ = newTargetGroups[g][i,3] * 1e6 - z_offset
            idx = findIdxOfClosestElem(adjustedZ, etlVals)
            newTargetGroups[g][i,3] = idx
        end
        if is1024
            newTargetGroups[g][:,1:2] *= 2
        else
            # @info "will alias fractional target center to nearest Int"
            # newTargetGroups[g] .= round.(newTargetGroups[g], digits=0)
        end
    end
    map(g -> all(isinteger.(g)) ? Int64.(g) : g, newTargetGroups)
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

function balanced_transition_order(nStims, nTransitionReps, max_attempts)
    best = Int64[]
    local successful
    for i=1:max_attempts
        trialOrder, successful = balanced_transition_order(nStims, nTransitionReps)
        if successful
            best = trialOrder
            break
        else
            if length(trialOrder) > length(best)
                println("new best: $(length(trialOrder))")
                best = trialOrder
            end
        end
    end
    if ~successful
        @warn "didn't succeed in making balanced trial order"
    end
    best, successful
end

"Try to find a trial order with balanced number of transitions."
function balanced_transition_order(nStims, nTransitionReps)
    nTrials = nStims^2 * nTransitionReps + 1

    transitionsLeft = ones(nStims,nStims) * nTransitionReps
    trialOrder = Int64[]
    successful = true
    push!(trialOrder, rand(1:nStims))
    for i in 2:nTrials
        prevStim = trialOrder[i-1]
        # @show prevStim
        # @show transitionsLeft
        candidateTransitionsLeft = findall(transitionsLeft[:,prevStim] .> 0)
        if length(candidateTransitionsLeft)==0
            successful=false
            break
        end
        nextStim = rand(candidateTransitionsLeft)
        # @show nextStim
        push!(trialOrder, nextStim)
        transitionsLeft[nextStim,prevStim] -= 1
    end
    trialOrder, successful
end

ez_gamma(x,gamma=0.5) = adjust_histogram(imadjustintensity(x), GammaCorrection(gamma))

# TODO: is there a better pattern to auto wrap..?
mutual_information(A,B) = py_utils.mutual_information(A,B)

# sadly, the below pattern will throw an error
# transparent_cmap = py_utils.transparent_cmap
# so we make a wrapper function instead...
# perhaps we should make a macro for this..?
function transparent_cmap(cmap; N=255, max_alpha=0.8)
    py_utils.transparent_cmap(cmap, N, max_alpha)
end