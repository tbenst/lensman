module DictThreadSafe
using Actors, Thunks
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

sv = @lazy supervisor()

# start dict server
function dictsrv(d::Dict; remote=false)
    dsrv = DictSrv(spawn(ds, d; remote))
    exec(dsrv.lk, Actors.supervise, reify(sv));
    dsrv
end

export DictSrv, dictsrv 

end