struct LazyTy5
    tyh5_path::String
    dset_str::String
    h5::HDF5.File
    dset::HDF5.Dataset
    function LazyTy5(tyh5_path, dset_str; mode="r", swmr=true)
        h5 = h5open(tyh5_path, mode, swmr=swmr)
        dset = h5[dset_str]
        new(tyh5_path, dset_str, h5, dset)
    end
end

"""
Index by HWZT & read from .ty.h5 file.

We ignore singleton Channel, and flip WH: WHCZT -> HWZT
"""
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

function Base.close(X::LazyTy5)
    close(X.h5)
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

abstract type AbstractDAG <: AbstractDict{Symbol, Any} end

struct DAG <: AbstractDAG
    nodes::Dict
    function DAG(;nodes...)
        new(nodes)
    end
end

function Base.getindex(r::AbstractDAG, k)
    reify(r.nodes[k])
end

function Base.setindex!(r::AbstractDAG, value, key...)
    r.nodes[key...] = value
end

Base.iterate(t::AbstractDAG) = iterate(t.nodes)
Base.iterate(t::AbstractDAG, i) = iterate(t.nodes, i)
Base.length(t::AbstractDAG) = length(t.nodes)