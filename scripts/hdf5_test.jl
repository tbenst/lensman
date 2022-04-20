using Pkg
Pkg.activate("..")
using Distributed, SharedArrays, ProgressMeter
rmprocs(workers()); addprocs(4); nworkers()
import Base.Threads: @threads
rand_h5_path = "/tmp/rand.h5"
# h5write(rand_h5_path, "/data", rand(100,100))
res = SharedArray{Float64}(10)
@everywhere begin
	using Pkg
	Pkg.activate("..")
	using HDF5, SharedArrays, ProgressMeter
end
# @showprogress @distributed for i in 1:10
@sync @distributed for i in 1:10
    h5open(rand_h5_path, "r", swmr=true) do thread_h5
        res[i] = sum(read(thread_h5["/data"]))
    end
end
println(sum(res))
