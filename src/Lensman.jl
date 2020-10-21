module Lensman

using AxisArrays, ANTsRegistration, NIfTI
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


export read_microns_per_pixel,
    read_mask,
    zbrain_units,
    antsApplyTransforms,
    getMaskNameIncludes
    # , segment_nuclei

end