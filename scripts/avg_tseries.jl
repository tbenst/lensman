using Lensman, FileIO, StatsBase

if ~size(ARGS,1)==1
    println(raw"Usage: julia --project=. scripts/avg_tseries.jl tifdir")
else
    tifDir = ARGS[1]
    tseries = loadTseries(tifDir);
    H, W, Z, T = size(tseries)
   save(tifDir*"_time-avg.tif", UInt16.(round.(mean(tseries, dims=4)[:,:,:,1], digits=0)))
end
