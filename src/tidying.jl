function per_trial_regions_df(tyh5_path, tseries_dset, window_len, stim_start_idx, stim_end_idx, trial_order,
    masks, region_names, stim_names=nothing
)
    nStimuli = maximum(trial_order)
    nTrials = length(stim_start_idx)
    if isnothing(stim_names)
        stim_names = string.(collect(1:nStimuli))
    end

    df = @sync @distributed vcat for i in 1:nTrials
        start_idx, end_idx = stim_start_idx[i], stim_end_idx[i]
        stim_name = stim_names[trial_order[i]]
        tempdf = DataFrame()
        _proc_region_df_f(tyh5_path, tseries_dset, masks, region_names, stim_name,
            start_idx, end_idx, i, window_len)
    end
    df
end

function _proc_region_df_f(tyh5_path, tseries_dset, masks, region_names, stim_name,
    start_idx, end_idx, trial, window_len, ϵ=0.0
)
    proc_tseries = LazyTy5(tyh5_path, tseries_dset)
    s = start_idx
    e = end_idx
    f0 = mean(proc_tseries[:,:,:,s-window_len:s-1],dims=4)[:,:,:,1]
    f = mean(proc_tseries[:,:,:,e+1:e+window_len],dims=4)[:,:,:,1]
    df = DataFrame()
    for (region_name, mask) in zip(region_names, masks)
        region_f = mean(f[mask])
        region_f0 = mean(f0[mask])
        region_response = (region_f - region_f0) / (region_f0 + ϵ)
        N = length(region_response)
        df = vcat(df, DataFrame("Δf/f"=>region_response,
            "region"=>region_name,
            "stim"=>stim_name,
            "trial"=>trial))
    end
    df
end