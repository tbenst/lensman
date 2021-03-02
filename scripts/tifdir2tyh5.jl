using Lensman

if ~size(ARGS,1)==2
    println(raw"Usage: julia --project=. scripts/tifdir2tyh5.jl tifdir z_offset")
    # println(raw"Usage: julia O:\users\tyler\code\lensman\notebooks\tifdir2tyh5.jl tifdir")
else
    tifDir = ARGS[1]
    z_offset = parse(Int, ARGS[2])
    slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
    stim_masks, stim_used_at_each_timestep = get_slm_stim_masks(tifDir, slmDir, z_offset)
    tseries = loadTseries(tifDir);
    write_experiment_to_tyh5(tseries, stim_masks, stim_used_at_each_timestep, tifDir*".ty.h5", compression_level=0)
end