using Lensman.DictThreadSafe
# using Dagger
# import Dagger: @par
# import Lensman.DictThreadSafe: dictsrv

# Recording(uri, settings) = Recording(uri, settings, dictsrv(Dict()))


# TODO: how do we handle the several formas of a tseries..?
# - h5 file dset
# - array
# string of h5_path for multiprocessing
# answer => use multiple dispatch + interface of required methods
# and LazyTy5 should be able to give strings for dset/h5path for multiprocessing


# TODO: split into functions that use a macro to fill args with same name..?
# eg suppose for this function:
# calc_trial_average(tseries::LazyTy5, stim_start_idx,
        # stim_end_idx, tseriesH, tseriesW, tseriesZ, trialOrder;
        # pre=16, post=16)
# we do:
# trial_average = @dagcall r calc_trial_average(tseries, stim_start_idx,
#        stim_end_idx, tseriesH, tseriesW, tseriesZ, trialOrder;
#        pre=16, post=16)
# which expands to:
# trial_average = calc_trial_average(r[:tseries], r[:stim_start_idx],
#        r[:stim_end_idx], r[:tseriesH], r[:tseriesW], r[:tseriesZ],
#        r[:trialOrder]; pre=16, post=16)
# ...but only if local variable isn't defined..? eg using @isdefined..?
# also would require new macro to replace @lazy such that:
# @lazy tseries_dir = find_folder(uri, tseries_root_dirs)
# becomes:
# tseries_dir = thunk(find_folder)(uri, tseries_root_dirs)
# r[:tseries_dir] = tseries_dir
# is this too much metaprogramming...? Harder for others to read...


"Master DAG for all computations on a recording."
function update_recording_dag(recording::DAG)
    @pun (uri, rel_plot_dir, tseries_dset, zseries_name, tseries_root_dirs,
        slm_dir, slm_root_dirs, lazy_tyh5, window_secs, zbrain_dir,
        oir_dir, warp_prefix, oir_920_name, oir_820_name, tyh5_path,
        h2b_zbrain, zbrain_units, rostral, dorsal, suite2p_dir, zbrain_masks,
        zbrain_mask_names
    ) = recording
    # procutil=Dict(Dagger.ThreadProc => 36.0)
    
    # TODO: figure out how to read slmTxtFile when multiple options...
    # (only needed for power fraction)
    # TODO: error doesn't show proper line number..
    # @par cache=true procutil=procutil begin
    # TODO: write new macro that first check if left hand symbol in
    # dict, if not then assign thunk. This way can rerun this function
    # multiple times to update DAG for live coding w/o invalidating cache
    @reversible begin
        tseries_dir = find_folder(uri, tseries_root_dirs)
        fish_dir = get_fish_dir(tseries_dir)
        zseries_dir = joinpath(fish_dir, zseries_name)
        expName = (x->splitpath(x)[end])(tseries_dir)
        recording_folder = (x->splitpath(x)[end-2])(tseries_dir)
        fish_name = (x->splitpath(x)[end-1])(tseries_dir)
        plot_dir = make_joinpath(fish_dir, rel_plot_dir)
        local_slm_dir = joinpath(fish_dir, "slm")
        tyh5_path = get_tyh5_path(tyh5_path, tseries_dir)
        test = count_threads()
        tseries = get_tseries(tseries_dir, tyh5_path, tseries_dset, lazy_tyh5)
        zseries = load_zseries(zseries_dir)
        zseries_xml = read_rec_xml(zseries_dir)
        tseries_xml = read_rec_xml(tseries_dir)
        zseries_zaxes = read_all_zaxis(zseries_xml)
        tseries_zaxis = read_first_zaxis(tseries_xml)   
        result = getExpData(tseries_xml)
        exp_date = ((x)->x[1])(result)
        frame_rate = ((x)->x[2])(result)
        etl_vals = ((x)->x[3])(result)
        res0 = tseriesTiffDirMetadata(tseries_dir)
        tseriesH = getindex(res0,1)
        tseriesW = getindex(res0,2)
        tseriesZ = getindex(res0,3)
        tseriesT = getindex(res0,4)        
        zseriesH = size(zseries,1)
        zseriesW = size(zseries,2)
        zseriesZ = size(zseries,3)
        vol_rate = frame_rate / tseriesZ
        date_str = Dates.format(exp_date, "dd-u-Y")
        tmp = find_folder(date_str, slm_root_dirs)
        used_slm_dir = (x->joinpath(splitpath(x)[1:end-1]...))(tmp)
        res = getTrialOrder(used_slm_dir, exp_date)
        trial_order = getindex(res,1)
        slm_exp_dir = getindex(res,2)
        first_target_group = ((x)->matread.(findMatGroups(x)[1]))(slm_exp_dir)
        slm_num = getSLMnum(first_target_group) # we assume same SLM for all groups
        power_per_cell = ((x)->first_target_group["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"])(first_target_group)
        res2 = slmpower(exp_date) # TODO: auto read from Oak..?
        slm1_power = getindex(res2,1)
        slm2_power = getindex(res2,2)
        voltageFile = glob_one_file("*VoltageRecording*.csv", tseries_dir)
        res3 = getStimTimesFromVoltages(voltageFile, tseriesZ)
        stim_start_idx = getindex(res3,1)
        stim_end_idx = getindex(res3,2)
        frameStartIdx = getindex(res3,3)
        nstim_pulses = getindex(res3,4)
        res4 = get_target_groups(slm_exp_dir)
        target_groups = getindex(res4,1)
        group_stim_freq = getindex(res4,1)
        nStimuli = maximum(trial_order)
        nTrials = size(trial_order,1)
        zOffset = getzoffset(exp_date, slm_num)
        suite2p_dir = isnothing(suite2p_dir) ? get_suite2p_dir(tseries_dir) : suite2p_dir
        nwb_path = joinpath(suite2p_dir, "ophys.nwb")
        ndict = read_nwb_rois(nwb_path)
        cell_traces = getindex(ndict,:cell_traces)
        cell_masks = getindex(ndict,:cell_masks)
        iscell = getindex(ndict,:iscell)
        nCells = getindex(ndict,:nCells)
        d0 = splitpath(zseries_dir)
        zseries_xml_path = (d->joinpath(d..., d[lastindex(d)] * ".xml"))(d0)
        zseries_xml = read_xml(zseries_xml_path);
        zseries_zaxes = read_all_zaxis(zseries_xml)
        # this is ugly, should fix with macro...
        et = ((e,t)-> e .+t)(etl_vals, tseries_zaxis)
        f = (x->z->searchsortedfirst(x, z))(zseries_zaxes)
        imaging2zseries_plane = map(f, et)
        window_len = ((vol_rate)->Int(round(window_secs * vol_rate)))(vol_rate)
        df_f_per_voxel_per_trial = get_df_f_per_voxel_per_trial_from_h5(
            tyh5_path, tseries_dset, stim_start_idx, stim_end_idx, window_len,
            tseriesH, tseriesW, tseriesZ)
        df_f_per_trial_dataframe = get_df_f_per_trial_dataframe(
            df_f_per_voxel_per_trial, trial_order)
        tseries_is1024 = (tseriesW==1024)
        twp1 = mapTargetGroupsToPlane(target_groups, etl_vals,
            is1024=tseries_is1024, zOffset=zOffset)
        targets_with_plane_index = map(x->Int.(round.(x, digits=0)), twp1)
        cells = makeCellsDF(targets_with_plane_index, stim_start_idx, stim_end_idx,
            trial_order, group_stim_freq)
        trial_average = calc_trial_average(tseries, stim_start_idx,
            stim_end_idx, tseriesH, tseriesW, tseriesZ, trial_order;
            pre=window_len, post=window_len)
        lateral_unit = microscope_lateral_unit(tseriesZ)
        target_size_px = spiral_size(exp_date, lateral_unit)
        zbrain_warps = glob("*SyN_Warped.nii.gz", fish_dir)
        zbrain_ants_cmd = (zbrain_warps -> length(zbrain_warps)==0
            ? ants_register(zseries, h2b_zbrain;
            interpolation = "WelchWindowedSinc", histmatch = 0,
            sampling_frac = 0.25, maxiter = 200, threshold=1e-8,
            use_syn = true, synThreshold = 1e-7, synMaxIter = 200,
            save_dir=fish_dir, dont_run = true)
            : "using cached ANTs for zbrain registration")(zbrain_warps)
        # TODO: should first try reading settings
        zbrain_warpedname = glob_one_file("*SyN_Warped.nii.gz", fish_dir)
        zbrain_transform_affine = glob_one_file("$warp_prefix*GenericAffine.mat", fish_dir)
        zbrain_transform_SyN = glob_one_file("$warp_prefix*SyN_1Warp.nii.gz", fish_dir)
        zbrain_transforms = [zbrain_transform_affine, zbrain_transform_SyN]
        zbrain_registered = niread(zbrain_warpedname)
        oir_920_file = joinpath(oir_dir, oir_920_name)
        oir_820_file = joinpath(oir_dir, oir_820_name)
        multimap_920 = read_olympus(oir_920_file)
        multimap_820 = read_olympus(oir_820_file)
        multimap_warps = glob("*SyN_Warped.nii.gz", fish_dir)
        multimap_ants_cmd = (multimap_warps -> length(multimap_warps)==0
            ? ants_register(zseries, multimap_920; interpolation = "WelchWindowedSinc",
            histmatch = 0, sampling_frac = 0.25, maxiter = 200, threshold=1e-8,
            use_syn = true, synThreshold = 1e-7, synMaxIter = 200,
            save_dir=fish_dir, dont_run = true) : "using cached ANTs for multimap registration")(zbrain_warps)
        multimap_warpedname = glob_one_file("$warp_prefix*SyN_Warped.nii.gz", oir_dir)
        mm920_registered = niread(multimap_warpedname)
        mm_transform_affine = glob_one_file("$warp_prefix*GenericAffine.mat", oir_dir)
        mm_transform_SyN = glob_one_file("$warp_prefix*SyN_1Warp.nii.gz", oir_dir)
        mm_transforms = [mm_transform_affine, mm_transform_SyN]
        mm820_registered = multimap_transforms(zseries, multimap_820,
            mm_transform_affine, mm_transform_SyN)
        region_mask_path = joinpath(fish_dir, "region_masks.h5")
        region_mask_dir = joinpath(fish_dir, "region_masks")
        region_masks_h5 = if isfile(region_mask_path)
            h5open(region_mask_path, "r", swmr=true)
        else
            save_region_masks(region_mask_path, zseries, zbrain_masks,
                zbrain_transforms;
                rostral=rostral, dorsal=dorsal)
        end

        sdict = read_suite2p(suite2p_dir)
        nCells = getindex(sdict,:nCells)
        cell_centers = getindex(sdict,:cell_centers)
        cells_mask = getindex(sdict,:cells_mask)
        iscell = getindex(sdict,:iscell)

        # targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etl_vals; is1024=(tseriesW==1024), zOffset=zOffset)
        # assigned all above

    end
    # if slm_num == 1
    #     slmpower_per_cell = @lazy slm1Power * power_per_cell / 1000
    # elseif slm_num == 2
    #     slmpower_per_cell = @lazy slm2Power * power_per_cell / 1000
    # end



    # procutil=Dict(Dagger.ThreadProc => 36.0)
    # @par cache=true procutil=procutil tseries =
    #     get_tseries(tseries_dir, tyh5_path, tseries_dset)

    @assign recording = (
        tseries_dir, fish_dir, zseries_dir, expName, recording_folder,
        fish_name, plot_dir, local_slm_dir, tyh5_path, test, tseries, zseries,
        zseries_xml, tseries_xml, zseries_zaxes, tseries_zaxis, result,
        exp_date, frame_rate, etl_vals, tseriesH, tseriesW, tseriesZ, tseriesT,
        zseriesH, zseriesW, zseriesZ, res, trial_order, slm_exp_dir,
        first_target_group, slm_num, power_per_cell, res2, slm1_power,
        slm2_power, voltageFile, res3, stim_start_idx, stim_end_idx, frameStartIdx,
        res4, target_groups, group_stim_freq, nStimuli, nTrials, zOffset,
        suite2p_dir, nwb_path, ndict, cell_traces, cell_masks, iscell, nCells,
        used_slm_dir, vol_rate, imaging2zseries_plane, zseries_xml,
        df_f_per_voxel_per_trial, df_f_per_trial_dataframe, lateral_unit,
        targets_with_plane_index, cells, trial_average, window_len,
        target_size_px, zbrain_registered, zbrain_ants_cmd,
        zbrain_registered, mm920_registered, mm820_registered, nstim_pulses,
        mm_transform_affine, mm_transform_SyN, zbrain_transforms,
        region_mask_path, zbrain_masks, region_masks_h5, zbrain_mask_names,
        nCells, cell_centers, cells_mask, iscell
    )

    recording
end

function update_dag(r)
    @pun (
        tseries_dir, fish_dir, zseries_dir, expName, recording_folder,
        fish_name, plot_dir, slm_dir, tyh5_path, test, tseries, zseries,
        zseries_xml, tseries_xml, zseries_zaxes, tseries_zaxis, result,
        exp_date, frame_rate, etl_vals, tseriesH, tseriesW, tseriesZ, tseriesT,
        zseriesH, zseriesW, zseriesZ, res, trial_order, slm_exp_dir,
        first_target_group, slm_num, power_per_cell, res2, slm1_power,
        slm2_power, voltageFile, res3, stim_start_idx, stim_end_idx, frameStartIdx,
        res4, target_groups, group_stim_freq, nStimuli, nTrials, zOffset,
        suite2p_dir, nwb_path, ndict, cell_traces, cell_masks, iscell, nCells
    ) = r.nodes

    if ~in(:somekey, keys(r.nodes))
        somekey = @lazy get_somekey()
        @assign r.nodes = somekey
    end

end

function get_target_groups(slmExpDir)
    target_groups = []
    group_stim_freq = []
    for mat in matread.(findMatGroups(slmExpDir))
        push!(target_groups, mat["cfg"]["maskS"]["targets"][1])
        push!(group_stim_freq, getMatStimFreq(mat))
    end
    target_groups, group_stim_freq
end

function count_threads()
    n = Threads.nthreads()
    a = zeros(n)
    Threads.@threads for i = 1:n
           a[i] = Threads.threadid()
    end
    length(unique(a))
end


function slow_job()
    sleep(5)
    10
end

"""Gather a Tuple/array of Dagger.Thunk

TODO: not very efficient as sequential, not parallel..?
"""
function gather(arr)
    # @par t = Tuple(arr...)
    # collect(t)
    map(collect, arr)
end

#### main analysis functions

"I thunk nothing."
ithunk() = thunk(identity)(nothing)

function is_node(nodes, key, continuation_thunk=ithunk())
    if key in keys(nodes)
        nodes[key]
    else
        continuation_thunk
    end
end

"Find a folder in first root_dirs where path exists."
function find_folder(name, root_dirs)
    for root_dir in root_dirs
        path = joinpath(root_dir, name)
        try
            if isdir(path)
                return path
            end
        catch
        end
    end
    @error "could not find folder $name"
end


function get_fish_dir(tseries_dir)
    sp = splitpath(tseries_dir)[1:end-1]
    joinpath(sp...)
end

function get_suite2p_dir(tseries_dir)
    sp = splitpath(tseries_dir)[1:end-1]
    joinpath(sp...)
end

function get_tyh5_path(settings_tyh5_path, tseries_dir)
    if isnothing(settings_tyh5_path)
        tseries_dir*".ty.h5"
    else
        settings_tyh5_path
    end
end

function get_tseries(tseries_dir, tyh5_path, tseries_dset, lazy_tyh5)
    if ~isnothing(tseries_dset)
        if isnothing(tyh5_path)
            tyh5_path = tseries_dir*".ty.h5"
        end
        if lazy_tyh5
            tseries = LazyTy5(tyh5_path, tseries_dset);
        else
            tseries = read_tyh5(tyh5_path, tseries_dset)
        end
    else
        tseries = loadTseries(tseries_dir);
    end
    tseries
end

function read_rec_xml(recording_dir)
    datafolders = splitpath(recording_dir)
    xmlPath = joinpath(datafolders..., datafolders[end] * ".xml")
    read_xml(xmlPath)
end