using SparseArrays, H5Sparse, Lensman

@testset "hdf5 threads" begin
    h5 = Lensman.HDF5ThreadSafe.h5server(tempname(), "w")
    s = sparse(randn(100,100) .> 0.95)
    s_vec = sparse(reshape(s,10000,1))
    H5SparseMatrixCSC(h5, "sparse", s_vec)
    s2 = reshape(H5SparseMatrixCSC(h5, "sparse"), 100,100)
    @test all(s2 .== s)
end

@testset "dictthreads" begin
    d = Lensman.DictThreadSafe.dictsrv(Dict{Int,Int}());
    d[1] = 0
    try
        d[:a] = 1 
    catch
    end
    d[2] = 1
    # @test d[2] == 1
    @assert d[2] == 1
    @assert d[1] == 0
end