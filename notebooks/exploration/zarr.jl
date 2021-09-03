using Zarr, HDF5, BenchmarkTools, Distributed
import Base.Threads: @threads

addprocs(36)
@everywhere using HDF5, Zarr
##
path = joinpath(mktempdir(), "test.zarr")
M, N = 1000, 10000
z1 = zcreate(Float32, M,N,path = path,chunks=(M/10, N/10))
# z1 = zcreate(Float32, M,N,path = path,chunks=(M/10, N/10))
m1 = rand(Float32, M,N);
z1[:,:] .= m1
h5_path = joinpath(mktempdir(), "test.h5")
h5 = h5open(h5_path, "w")
h5["a"] = m1;
close(h5)
h5 = h5open(h5_path, "r"; swmr=true)

##

function threads_sum(arr)
    result = zeros(size(arr,1))

    @threads for i = 1:size(arr,1)
        result[i] = sum(arr[i,:])
    end
    result
end


@everywhere function distributed_sum(arr::HDF5.File; dset="a")
    @distributed vcat for i = 1:size(arr[dset],1)
        h5 = h5open(arr.filename)
        array = h5[dset]
        sum(array[i,:])
        close(h5)
    end
end

@everywhere function distributed_sum(arr::ZArray)
    @distributed vcat for i = 1:size(arr,1)
        sum(arr[i,:])
    end
end

function seq_sum(arr)
    result = zeros(size(arr,1))

    for i = 1:size(arr,1)
        result[i] = sum(arr[i,:])
    end
    result
end
##
@btime seq_sum(m1); # 44ms
@btime threads_sum(m1); # 3ms
##
@btime seq_sum(h5["a"]); # 4.37s
@btime seq_sum(z1); # 2.08s
@btime threads_sum(z1); # 3.36s
@btime distributed_sum(h5) # 194ms
@btime distributed_sum(z1) # 590ms