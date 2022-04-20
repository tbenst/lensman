if ~size(ARGS,1)==1
    println(raw"Usage: julia --project=. scripts/artifact_filter.jl uri [outname]")
    exit()
end

uri = ARGS[1]
if length(ARGS) > 1
    outname = ARGS[2]
    global use_outname = true
else
    global use_outname = false
end
using Lensman
import Lensman: @pun
resources = Resources();
recording = Recordings[uri](;resources...,
    # channel_str="Ch2"
    tseries_read_strategy=:lazy_tiff
);

@pun (fish_dir, exp_name, tseries, tseriesT, artifacts
) = recording;

# kalman_path = joinpath(fish_dir, exp_name * "_red_kalman.h5")
if use_outname
    outpath = outname
else
    outpath = joinpath(fish_dir, exp_name * "_no_artifact.h5")
    outpath = replace(outpath , "/scratch" => "/data/dlab")
end
mkpath(dirname(outpath))
Lensman.artifact_filter_stream(outpath, tseries, tseriesT, artifacts);
println(outpath)