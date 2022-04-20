using HDF5, H5Sparse, SparseArrays
import Base.Threads: @threads
function write_n_sparse_datasets(h5path, n, channel, blocker, type=:H5)
    h5 = h5open(h5path, "w")
    while n > 0
        @assert Threads.threadid() == 1 # we only access from the master thread
        dset_name, data = take!(channel)
        take!(blocker)
        println("take $n")
        if type == :H5
            h5[dset_name] = data
        elseif type == :H5Sparse
            H5SparseMatrixCSC(h5, dset_name, data)
        end
        n -= 1
    end
    h5
end

function test(type=:H5, n=10, N=512*512*300)
    to_write = Channel(Inf)
    blocker = Channel(5) # restrict to N active threads
    h5_path = tempname()*".h5"
    isfile(h5_path) ? rm(h5_path) : ()
    @threads for i in 1:n
        @async begin
            put!(blocker, i)
            println("start $i")
            data = collect(rand(N,1) .> 0.9)
            if type == :H5
                put!(to_write, "$i" => data)
            elseif type == :H5Sparse
                data = sparse(data)
                put!(to_write, "$i" => data)
            end
            println("put $i")
        end
    end
    println("call write")
    h5 = write_n_sparse_datasets(h5_path, n, to_write, blocker, type)

    @show keys(h5)
    close(h5)
    rm(h5_path)
end

test(:H5Sparse, 50) # always errors