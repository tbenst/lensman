using Lensman, HDF5
if ~size(ARGS,1)==2
    println(raw"Usage: julia --project=. scripts/tyh52suite2ph5.jl .ty.h5")
else
    tyh5 = ARGS[1]
    tseries = loadTseries(tyh5);
    outname = tyh5*"_s2p.h5"
    h5 = h5open(outname, "cw")
    h5["/data"] = tseries;
    close(h5)
end
