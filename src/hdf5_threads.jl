# thread-safe write access to HDF5 using Actors
module HDF5ThreadSafe

using Actors, HDF5
import Actors: spawn

struct HDF5server{L}
    lk::L
end
(h5server::HDF5server)(f::Function, args...) = call(h5server.lk, f, args...)
(h5server::HDF5server)() = call(h5server.lk)

# indexing interface
Base.getindex(h5s::HDF5server, key) = call(h5s.lk, getindex, key)
Base.setindex!(h5s::HDF5server, value, key) = call(h5s.lk, setindex!, value, key)
Base.close(h5s::HDF5server) = call(h5s.lk, close)

# hdf5 server behavior
h5server(h5::HDF5.File, f::Function, args...) = f(h5, args...)
h5server(h5::HDF5.File) = show(h5)

# start hdf5 server (constructor)
function hdf5server(filename::AbstractString, mode::AbstractString="r"; swmr=false,
    remote=false)
    h5 = h5open(filename, mode; swmr=swmr)
    HDF5server(spawn(h5server, h5; remote))
end

export HDF5server, hdf5server 

end