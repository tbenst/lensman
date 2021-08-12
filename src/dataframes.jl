"df/f per trial per voxel tidy dataframe."
function get_df_f_per_trial_dataframe(df_f_per_voxel_per_trial, trial_order)
    nTrials = size(df_f_per_voxel_per_trial,1)
    nVoxels = prod(size(df_f_per_voxel_per_trial)[2:end])
    df = DataFrame(reshape(df_f_per_voxel_per_trial, nTrials, nVoxels),
            map(string,collect(1:nVoxels)))
    df[!, :stim] = trial_order
    df = stack(df, 1:nVoxels)
    rename!(df, Dict(:variable => :voxel_id, :value => :df_f))
    df[!,:voxel_id] = parse.(Int,df[!,:voxel_id])
    df
end