# thread-safe write access to HDF5 using Actors
module HDF5ThreadSafe

using Actors, HDF5, Thunks
import Actors: spawn
import HDF5: create_group, create_dataset

struct HDF5server{L} <: HDF5.H5DataStore
    id::HDF5.hid_t
    filename::String
    lk::L
end
(_h5server::HDF5server)(f::Function, args...) = call(_h5server.lk, f, args...)
(_h5server::HDF5server)() = call(_h5server.lk)

# indexing interface
Base.getindex(h5s::HDF5server, key) = call(h5s.lk, getindex, key)
Base.setindex!(h5s::HDF5server, value, key) = call(h5s.lk, setindex!, value, key)
Base.close(h5s::HDF5server) = call(h5s.lk, close)
Base.keys(h5s::HDF5server) = call(h5s.lk, keys)
create_group(h5s::HDF5server, string) = call(h5s.lk, create_group, string)
create_dataset(h5s::HDF5server, args; kwargs...) = call(h5s.lk, create_dataset, args...; kwargs...)
delete_object(h5s::HDF5server, args; kwargs...) = call(h5s.lk, delete_object, args...; kwargs...)

# hdf5 server behavior
_h5server(h5::HDF5.File, f::Function, args...) = f(h5, args...)
_h5server(h5::HDF5.File) = call(h5s.lk, show)

sv = @lazy supervisor()

# start hdf5 server (constructor)
"""Open or create an HDF5 file using an Actors server for parallel read/write.

See HDF5.h5open documentation for more information on modes.
"""
function h5server(filename::AbstractString, mode::AbstractString="r"; swmr=true,
    remote=false)
    h5 = h5open(filename, mode; swmr=swmr)
    fn = h5.filename
    id = h5.id
    hsrv = HDF5server(id, fn, spawn(_h5server, h5; remote))
    exec(hsrv.lk, Actors.supervise, reify(sv))
    hsrv
end

"""Support do-block syntax.

Based on:
https://docs.julialang.org/en/v1/manual/functions/#Do-Block-Syntax-for-Function-Arguments
"""
function h5server(f::Function, args...)
    h5 = h5server(args...)
    try
        f(h5)
    finally
        close(h5)
    end
end


export HDF5server, h5server

end