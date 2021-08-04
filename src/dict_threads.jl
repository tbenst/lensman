module DictThreadSafe
using Actors
import Actors: spawn

struct DictSrv{L}
    lk::L
end
(ds::DictSrv)(f::Function, args...) = call(ds.lk, f, args...)
(ds::DictSrv)() = call(ds.lk)

# indexing interface
Base.getindex(d::DictSrv, key) = call(d.lk, getindex, key)
Base.setindex!(d::DictSrv, value, key) = call(d.lk, setindex!, value, key)
Base.keys(d::DictSrv) = call(d.lk, keys)

# dict server behavior
ds(d::Dict, f::Function, args...) = f(d, args...)
ds(d::Dict) = copy(d)

# start dict server
dictsrv(d::Dict; remote=false) = DictSrv(spawn(ds, d; remote))

export DictSrv, dictsrv 

end