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


function per_trial_regions_df(tseries::LazyTy5, window_len, stim_start_idx, stim_end_idx, trial_order,
    masks, region_names, stim_names=nothing
)
    per_trial_regions_df(tseries.tyh5_path, tseries.dset_str, window_len, stim_start_idx, stim_end_idx,
        trial_order, masks, region_names, stim_names)
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

"Make dataframe showing first f1_end seconds."
function tseries_ramping_df(tseries, genotype, fish; vol_rate = 3, f0_end=6,
     f1_end=72, ϵ=1
)
    df = DataFrame()
    f0 = mean(tseries[:,:,:,1:f0_end])
    f = mean(tseries[:,:,:,1:f1_end], dims=[1,2,3])[1,1,1,:]
    rel_time = collect(1:length(f))
    time = rel_time./vol_rate
    df_f = @. (f-f0)/(f0 .+ ϵ)
    df = DataFrame(
        rel_time=rel_time,
        time= time,
        df_f= df_f,
        f = f,
        genotype=genotype,
        fish=fish
    )
end

"Make dataframe showing first f1_end seconds."
function recording_to_df_df(recording, genotype)
    try
        @pun (tseries, vol_rate, uri) = recording
        # f0_end = Int(round(vol_rate * 2))
        f0_end=1
        f1_end = Int(round(vol_rate * 20))-1
        fish = uri
        tseries_ramping_df(tseries, genotype, fish; vol_rate = vol_rate, f0_end=f0_end,
            f1_end=f1_end)
    catch
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println(stdout)
        end
        DataFrame()
    end
end


function add_period_to_df(df)
    df[!,:period] = ["" for _ in 1:size(df,1)];
    @assert maximum(df[!,:stim]) == 16
    for s in 1:16
        if s <= 6
            p = "early"
        elseif s <= 11
            p = "mid"
        else
            p = "late"
        end
        idxs = df[!, :stim] .== s
        df[idxs,:period] .= p
    end
    catarr = CategoricalArray(df[!,:period], ordered=true)
    levels!(catarr, ["early", "mid", "late"])
    df[!,:period] = catarr
    df
end

function add_major_regions_to_df(df)
    df[!,"major region"] = ["" for _ in 1:size(df,1)];
    rhom_regions = unique(filter(r->occursin("Rhombencephalon", r), df.region))
    max = 15
    times, remainder = divrem(length(rhom_regions),max)
    nrhom_sub = times + (remainder > 0)
    rhom_map = Dict(r => "Rhombencephalon $(Char((i % nrhom_sub) + 65))"
        for (i,r) in enumerate(rhom_regions))
    for i in 1:size(df,1)
        # @show df[i,:region]
        if df[i,:region] == "Spinal Cord"
            df[i,"major region"] = "Spinal Cord"
        else
            m = match(r"^([\w\s]+) -\s?(.*)", df[i,:region])
            if ~isnothing(m)
                if m[1] == "Rhombencephalon"
                    # too many, split into subregions
                    df[i,"major region"] = rhom_map[df[i,:region]]
                else
                    df[i,"major region"] = m[1]
                end
            else
                error("no match for $(df[i,:region])")
            end
        end
    end
    df
end

function rename_regions_in_df(df)
    for i in 1:size(df,1)
        # @show df[i,:region]
        if df[i,:region] == "Spinal Cord"
            df[i,"major region"] = "Spinal Cord"
        else
            m = match(r"^([\w\s]+) -\s?(.*)", df[i,:region])
            if ~isnothing(m)
                df[i,:region] = m[2]
            else
                error("no match for $(df[i,:region])")
            end
        end
    end
    df
end
