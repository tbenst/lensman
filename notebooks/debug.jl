using HDF5, H5Sparse, SparseArrays
import Base.Threads: @threads
function write_n_datasets(h5path, n, channel,blocker)
    h5 = h5open(h5path, "w")
    while n > 0
        dset_name, data = take!(channel)
        take!(blocker)
        println("take $n")
        H5SparseMatrixCSC(h5, dset_name, data)
        n -= 1
    end
    close(h5)
end

n = 10
to_write = Channel(Inf)
blocker = Channel(5) # restrict to N active threads
h5_path = tempname()*".h5"
isfile(h5_path) ? rm(h5_path) : ()
@threads for i in 1:n
    @async begin
        put!(blocker, i)
        println("start $i")
        data = collect(rand(10,10) .> 0.9)
        data = sparse(data)
        sleep(1)
        put!(to_write, ("$i", data))
        println("put $i")
    end
end
println("call write")
# write_n_datasets(h5_path, n, to_write, nothing)
write_n_datasets(h5_path, n, to_write, blocker)


# h5open(h5_path, "r") do h5
#     @show h5
# end
h5 = h5open(h5_path, "r")
@show keys(h5)
close(h5)
##








##
using Actors, HDF5
import Actors: spawn, newLink
import Base.Threads: @threads

test()  = 1

_calc(addr::Link, f::Function, args, kwargs) = send(addr, f(args...;kwargs...))
_calc(addr::Link, f::Function, args) = send(addr, f(args...))
_calc(addr::Link, f::Function) = send(addr, f())
actor_responses = newLink()
actor = spawn(_calc, actor_responses)
sv = supervisor()
exec(actor_responses, supervise, sv)
"Blocking interface to the actor."
function calc(f, args...; kwargs...)
    # request(actor, f, args, kwargs; timeout=Inf)[1]
    send(actor, f, args, kwargs; timeout=Inf)[1]
end
##
x = calc(test)
##
h5 = calc(h5open, tempname()*".h5", "w")
##
function assign_idx(h5, i)
    h5["$i"] = rand(10,10)
end

@threads for i in 1:10
    calc(assign_idx, h5, i)
end
h5
##
