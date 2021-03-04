using Lensman

# note 
if ~size(ARGS,1)==1
    println(raw"Usage: julia --project=. scripts/tifdir2tyh5.jl tifdir")
    println("for windows, try doing `set JULIA_NUM_THREADS=12` first?")
    #  set JULIA_NUM_THREADS=12
    # println(raw"Usage: julia O:\users\tyler\code\lensman\notebooks\tifdir2tyh5.jl tifdir")
else
    tifDir = ARGS[1]
    tseries = loadTseries(tifDir);
    write_experiment_to_tyh5(tseries, tifDir*".ty.h5", compression_level=0)
end