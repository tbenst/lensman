using Lensman, HDF5
if ~size(ARGS,1)==2
    println(raw"Usage: julia --project=. scripts/bruker2suite2ph5.jl tifdir")
else
    tifDir = ARGS[1]
    tseries = loadTseries(tifDir);
    outname = tifDir*"_s2p.h5"
    h5 = h5open(outname, "w")
    h5["/data"] = tseries;
    close(h5)
end
