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


struct LazyTiff
    framePlane2tiffPath::Dict{Tuple, String}
    tseries_size::NTuple{4, Int}
    function LazyTiff(tiff_dir; contains_str="Ch3")
        H, W, Z, T, framePlane2tiffPath = tseriesTiffDirMetadata(tiff_dir, contains_str)
        new(framePlane2tiffPath, (H,W,Z,T))
    end
end

"""
Index by HWZT & read from .ty.h5 file.

We ignore singleton Channel, and flip WH: WHCZT -> HWZT
"""
function Base.getindex(X::LazyTy5,i...)
    h,w,z,t = i
    c = 1
    # res = Base.getindex(X.dset, w, h, z, c, t)
    res = X.dset[w, h, z, c, t]
    nd = ndims(res)
    if nd > 0
        permutedims(res, (2,1,collect(3:nd)...))
    else
        res
    end
end

# TODO: try ReusePatterns to make this DRY

function Base.size(X::LazyTy5)
    W, H, Z, C, T = size(X.dset)
    H, W, Z, T
end

# function Base.axes(X::LazyTiff, args...)
#     Base.axes(X.dest, args...)
# end

function Base.axes(X::LazyTy5, args...)
    Base.axes(X.dest, args...)
end


function Base.axes(X::Union{LazyTy5,LazyTiff},i)
    # Base.axes(Base.OneTo(size(X,i)))
    Base.axes(X.dset,i)
end

function Base.close(X::LazyTy5)
    close(X.h5)
end

function Base.size(X::LazyTy5, i)
    W, H, Z, C, T = size(X.dset)
    (H, W, Z, T)[i]
end

"""
"""
function Base.getindex(X::LazyTiff,i...)
    Hs,Ws,Zs,Ts = i
    H, W = size(X)[1:2]
    if typeof(Zs) == Colon
        Z = size(X,3)
        Zs = 1:Z
    else
        Z = length(Zs)
    end
    if typeof(Ts) == Colon
        T = size(X,4)
        Ts = 1:T
    else
        T = length(Ts)
    end

    # check if iterable
    if ~applicable(iterate, Ts)
        Ts = [Ts]
    end
    if ~applicable(iterate, Zs)
        Zs = [Zs]
    end

    tseries = Array{UInt16}(undef, H, W, Z, T)
    # @threads for (i,t) in enumerate(Ts)
    for (i,t) in enumerate(Ts)
        for (j,z) in enumerate(Zs)
            tp = X.framePlane2tiffPath[t,z]
            tseries[:,:,j,i] .= reinterpret(UInt16, ImageMagick.load(tp))
        end
    end
    tseries[Hs, Ws, :, :]
end


function Base.size(X::LazyTiff)
    X.tseries_size
end


function Base.size(X::LazyTiff, i)
    X.tseries_size[i]
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