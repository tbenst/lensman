if ~size(ARGS,1)==1
    println(raw"Usage: julia --project=. scripts/avg_tseries.jl uri")
    exit()
end

using Lensman
import Lensman: @pun
@wan "using red channel"
resources = Resources();
recording = Recordings[uri](;resources...,
    # channel_str="Ch2"
);

@pun (fish_dir, exp_name, tseries, tseriesT, artifacts
) = recording;

# kalman_path = joinpath(fish_dir, exp_name * "_red_kalman.h5")
kalman_path = joinpath(fish_dir, exp_name * "kalman.h5")
kalman_path = replace(kalman_path , "/scratch" => "/data/dlab")
mkpath(dirname(kalman_path))
Lensman.kalman_filter_stream(kalman_path, tseries, tseriesT, artifacts);
kalman_path