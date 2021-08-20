function per_trial_regions_df(tseries::Array, window_len, stim_start_idx, stim_end_idx, trial_order,
    masks, region_names, stim_names=nothing
)
    nStimuli = maximum(trial_order)
    nTrials = length(stim_start_idx)
    if isnothing(stim_names)
        stim_names = string.(collect(1:nStimuli))
    end

    df = vcat(collect(Folds.map(1:nTrials) do i
        start_idx, end_idx = stim_start_idx[i], stim_end_idx[i]
        stim_name = stim_names[trial_order[i]]
        _proc_region_df_f(tseries, masks, region_names, stim_name,
            start_idx, end_idx, i, window_len)
    end)...)
    df
end

function per_trial_regions_df(tyh5_path::AbstractString, tseries_dset::AbstractString, window_len, stim_start_idx, stim_end_idx, trial_order,
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
        _proc_region_df_f(tyh5_path, tseries_dset, masks, region_names, stim_name,
            start_idx, end_idx, i, window_len)
    end
    df
end

"Return DataFrame of df/f values"
function imap_DF_by_region(imap, region_masks_h5::HDF5.File, regions, imaging2zseries_plane,
    stim_names=nothing, pretty_region_names=nothing
)
    nStimuli = size(imap,4)
    if isnothing(stim_names)
        stim_names = string.(collect(1:nStimuli))
    end
    pretty_region_names = isnothing(pretty_region_names) ? regions : pretty_region_names

    H, W, Z = size(imap)[1:3]
        
    df = DataFrame()
    for (pretty, region) in zip(pretty_region_names, regions)
        mask = read_registered_mask(region_masks_h5, region)[:,:,imaging2zseries_plane]
        for s in 1:nStimuli
            voxel_responses = selectdim(imap,4,s)[mask]
            voxel_responses = mean(voxel_responses)
            N = length(voxel_responses)
            df = vcat(df, DataFrame("Δf/f"=>voxel_responses,
                # region=repeat(CategoricalArray([region]),N),
                "region"=>pretty,
                # stim=repeat(CategoricalArray([stim_names[s]]),N)))
                "stim"=>stim_names[s]))
        end
    end
    df
end

function imap_DF_by_region(imap, masks, regions, imaging2zseries_plane,
    H, W, Z, stim_names=nothing
)
    nStimuli = size(imap,4)
    if isnothing(stim_names)
        stim_names = string.(collect(1:nStimuli))
    end

    H, W, Z = size(imap)[1:3]
        
    df = DataFrame()
    for (region, mask) in zip(regions, masks)
        for s in 1:nStimuli
            voxel_responses = selectdim(imap,4,s)[mask]
            N = length(voxel_responses)
            df = vcat(df, DataFrame("Δf/f"=>voxel_responses,
                "region"=>region,
                "stim"=>stim_names[s]))
        end
    end
    df
end

function _proc_region_df_f(tyh5_path::AbstractString, tseries_dset::AbstractString, masks,
    region_names, stim_name, start_idx, end_idx, trial, window_len, ϵ=0.0
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
        if occursin("left ", region_name)
            hemisphere="left"
            region_name = replace(region_name, "left " => "")
        elseif occursin("right ", region_name)
            hemisphere="right"
            region_name = replace(region_name, "right " => "")
        else
            hemisphere="unknown"
        end
        df = vcat(df, DataFrame("Δf/f"=>region_response,
            "region"=>region_name,
            "stim"=>stim_name,
            "trial"=>trial,
            "hemisphere"=>hemisphere))
    end
    df
end

function _proc_region_df_f(tseries::Array, masks, region_names, stim_name,
    start_idx, end_idx, trial, window_len, ϵ=0.0
)
    s = start_idx
    e = end_idx
    f0 = mean(tseries[:,:,:,s-window_len:s-1],dims=4)[:,:,:,1]
    f = mean(tseries[:,:,:,e+1:e+window_len],dims=4)[:,:,:,1]
    df = DataFrame()
    for (region_name, mask) in zip(region_names, masks)
        region_f = mean(f[mask])
        region_f0 = mean(f0[mask])
        region_response = (region_f - region_f0) / (region_f0 + ϵ)
        N = length(region_response)
        if occursin("left", region_name)
            hemisphere="left"
            region_name = replace(region_name, "left " => "")
        elseif occursin("right", region_name)
            hemisphere="right"
            region_name = replace(region_name, "right " => "")
        else
            hemisphere="unknown"
        end
        df = vcat(df, DataFrame("Δf/f"=>region_response,
            "region"=>region_name,
            "stim"=>stim_name,
            "trial"=>trial,
            "hemisphere"=>hemisphere))
    end
    df
end