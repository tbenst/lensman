struct LazyTy5
    dset::HDF5.Dataset
end

"Index by HWZT & read from .ty.h5 file."
function Base.getindex(X::LazyTy5,i...)
    h,w,z,t = i
    c = 1
    res = Base.getindex(X.dset, w, h, z, c, t)
    nd = ndims(res)
    if nd > 0
        permutedims(res, (2,1,collect(3:nd)...))
    else
        res
    end
end

function Base.size(X::LazyTy5)
    W, H, Z, C, T = size(X.dset)
    H, W, Z, T
end

function Base.size(X::LazyTy5, i)
    W, H, Z, C, T = size(X.dset)
    (H, W, Z, T)[i]
end

"Take Float32 [0,1], and convert to UInt16"
function float2uint(x)
    # PV max is 8192..? or at least, that's what we used in babelfish
    ret = round.(x .* 8192, digits=0)
    clamp!(ret,0,2^16-1)
    convert(Array{UInt16},ret)
end
