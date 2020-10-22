module Lensman

using AxisArrays, ANTsRegistration, NIfTI, ImageMagick, Images,
    ImageDraw, ImageFiltering, PyCall, MAT, Dates, DataStructures
import Base.Threads.@threads

# https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
peak_local_max = PyNULL()
disk = PyNULL()
match_template = PyNULL()

function __init__()
    pyimport("skimage")
    copy!(peak_local_max, pyimport("skimage.feature").peak_local_max)
    copy!(disk, pyimport("skimage.morphology").disk)
    copy!(match_template, pyimport("skimage.feature").match_template)
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
function findNeurons(plane)
    rescaled_img = (adjust_histogram(imadjustintensity(plane), Equalization()))
    # mask out background
    laplace_mask = opening(imadjustintensity(morpholaplace(rescaled_img)).<0.5)
    Gray.(laplace_mask)
    ##
    img = rescaled_img .* laplace_mask
    feat_img = match_template(Float64.(img), disk(4), pad_input=true, mode="constant", constant_values=0)
    local_maxi = peak_local_max(feat_img, indices=false, footprint=disk(4), exclude_border=false, threshold_rel=.1)
end

"Make .mat file for Sean's multiSLM stim software."
# ::Array{Float64,2}
function create_targets_mat(targets, outname::String;
    outdir = "/mnt/deissero/users/tyler/slm/masks/")
    today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
    today_dir = mkpath(joinpath(outdir, today_str))
    matpath = joinpath(today_dir, "$outname.mat")

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

    out_mat["targets"] = Array{Any}(nothing,1,1)
    out_mat["targets"][1] = copy(targets) # ensure not an Adjoint type
    out_mat["zrange"] = [minimum(targets[:,3]) maximum(targets[:,3])]
    # channel 1 (red) and 2 (green), used for GUI preview

    ntargets = size(targets,1)
    # targetsGroup must alternate between Float32(1,3) and zeros(0,0), and end with 2 * zeros(0,0)
    out_mat["targetsGroup"] = Array{Any}(nothing,1,ntargets*2+2)
    for (i, v) in enumerate(eachrow(targets))
        # SubArrays cannot be written, so we convert to Array
        out_mat["targetsGroup"][i*2-1] = convert(Array,v')
        out_mat["targetsGroup"][i*2] = zeros(0,0)
    end
    out_mat["targetsGroup"][end-1] = zeros(0,0)
    out_mat["targetsGroup"][end] = zeros(0,0)

    matwrite(matpath, stim_struct)
    # strip /mnt/deissero...
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

function create_slm_stim(target_groups, outname::String)
    targets_mats = []
    for (i, targets) in enumerate(target_groups)
        name = "$(outname)_group_$i"
        push!(targets_mats,create_targets_mat(targets, name))
    end
    trials_txt = create_trials_txt(targets_mats, outname)
    return targets_mats, trials_txt
end

"Convert mask to an Array of indices, similar to argwhere in Python."
function mask2ind(mask)
    hcat(collect.(Tuple.(findall(mask)))...)'
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
    mask2ind
    # , segment_nuclei

end