using HDF5
import Base.Threads: @threads
function write_n_datasets(h5path, n, channel, blocker, type=:H5)
    h5 = h5open(h5path, "w")
    while n > 0
        @assert Threads.threadid() == 1 # we only access from the master thread
        dset_name, data = take!(channel)
        println("take $n")
        take!(blocker)
        if type == :H5
            h5[dset_name] = data
        elseif type == :create_group
            group = create_group(h5, dset_name)
            group[dset_name] = data
        elseif type == :create_dset
            group = create_group(h5, dset_name)
            create_dataset(
                group, "data", 
                eltype(data), 
                ((length(data),), (-1,)),
                chunk=HDF5.heuristic_chunk(data),
                blosc=5)[1:length(data)] = data
        elseif type == :simple_create_dset
            group = create_group(h5, dset_name)
            create_dataset(group, "data", data)
        elseif type == :copy_create_dset
            data = copy(data)
            dset_name = string(dset_name)
            group = create_group(h5, dset_name)
            create_dataset(group, "data", data)
        end
        n -= 1
    end
    h5
end

function test(type=:H5, n=10, N=512*512*300)
    to_write = Channel(Inf)
    blocker = Channel(5)
    h5_path = tempname()*".h5"
    isfile(h5_path) ? rm(h5_path) : ()
    @threads for i in 1:n
        @async begin
            put!(blocker, i)
            println("start $i")
            data = collect(rand(N,1) .> 0.9)
            if type in [:H5, :create_group]
                put!(to_write, "$i" => data)
            elseif type in [:create_dset, :simple_create_dset, :copy_create_dset]
                data = reshape(data, N)
                put!(to_write, "$i" => data)
            else
                error("bad type")
            end
            println("put $i")
        end
    end
    println("call write")
    h5 = write_n_datasets(h5_path, n, to_write, blocker, type)

    @show keys(h5)
    close(h5)
    rm(h5_path)
    println("succeeded")
end

# test(:H5, 10) # always works
# test(:H5, 50) # always works
# # test(:H5, 500) # takes a long time, but works
# test(:create_group, 100) # no error
# test(:create_dset, 10) # no error
# # test(:create_dset, 100) # usually errors
# test(:simple_create_dset, 5) # usually succeeds
# test(:simple_create_dset, 10) # sometimes errors
# test(:simple_create_dset, 100) # usually errors
test(:copy_create_dset, 100) # usually errors