using Lensman.DictThreadSafe
using Dagger
import Dagger: @par
using Thunks, NPZ
# import Lensman.DictThreadSafe: dictsrv
struct Recording
    uri::String
    settings::Dict
    thunks::Dict
    Recording(uri; settings...) = new(
        uri, settings, make_dag(uri, merge(DEFAULT_SETTINGS, settings)))
end


# Recording(uri, settings) = Recording(uri, settings, dictsrv(Dict()))
DEFAULT_SETTINGS = Dict(
    :tseries_root_dirs => [
        "/oak/stanford/groups/deissero/users/tyler/b115",
        "/scratch/b115",
        "/scratch2",
        "/data/dlab/b115",
        "/mnt/deissero/users/tyler/b115"
    ],
    :slm_root_dirs => [
        "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files",
        "/mnt/deissero/users/tyler/b115/SLM_files",
        "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/",
        "/mnt/deissero/users/tyler/b115/SLM_files/",
        "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/"
    ],
    :rel_plot_dir => "plots",
    :tseries_dset => nothing,
    :tyh5_path => nothing,
    :zbrain_units => (0.798μm, 0.798μm, 2μm),
    :slm_dir => "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/",
    :suite2p_dir => nothing,
    # slm_dir => "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files",
    # slm_dir => "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files",
    # slm_dir => "/mnt/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment",
)

"Master DAG for all computations."
function make_dag(uri, settings)
    dag = Dict()
    @pun (rel_plot_dir, tseries_dset, zseries_name, tseries_root_dirs,
        slm_dir, slm_root_dirs
    ) = settings
    # procutil=Dict(Dagger.ThreadProc => 36.0)
    
    # TODO: figure out how to read slmTxtFile when multiple options...
    # (only needed for power fraction)
    # TODO: error doesn't show proper line number..
    # @par cache=true procutil=procutil begin
    # TODO: write new macro that first check if left hand symbol in
    # dict, if not then assign thunk. This way can rerun this function
    # multiple times to update DAG for live coding w/o invalidating cache
    @thunk begin
        tseries_dir = find_folder(uri, tseries_root_dirs)
        fish_dir = get_fish_dir(tseries_dir)
        zseries_dir = joinpath(fish_dir, zseries_name)
        expName = (x->splitpath(x)[end])(tseries_dir)
        recording_folder = (x->splitpath(x)[end-2])(tseries_dir)
        fish_name = (x->splitpath(x)[end-1])(tseries_dir)
        plot_dir = joinpath(fish_dir, rel_plot_dir)
        local_slm_dir = joinpath(fish_dir, "slm")
        tyh5_path = get_tyh5_path(settings[:tyh5_path], tseries_dir)
        test = count_threads()
        tseries = get_tseries(tseries_dir, tyh5_path, tseries_dset)
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
        stimStartIdx = getindex(res3,1)
        stimEndIdx = getindex(res3,2)
        frameStartIdx = getindex(res3,3)
        res4 = get_target_groups(slm_exp_dir)
        target_groups = getindex(res4,1)
        group_stim_freq = getindex(res4,1)
        nStimuli = maximum(trial_order)
        nTrials = size(trial_order,1)
        zOffset = getzoffset(exp_date, slm_num)
        suite2p_dir = is_setting(settings, :suite2p_dir, thunk(get_suite2p_dir)(tseries_dir))
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
        # sdict = read_suite2p(suite2p_dir)
        # nCells = getindex(sdict,:nCells)
        # cell_centers = getindex(sdict,:cell_centers)
        # cells_mask = getindex(sdict,:cells_mask)

        # targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etl_vals; is1024=(tseriesW==1024), zOffset=zOffset)
        # assigned all above

    end

    # if slm_num == 1
    #     slmpower_per_cell = @thunk slm1Power * power_per_cell / 1000
    # elseif slm_num == 2
    #     slmpower_per_cell = @thunk slm2Power * power_per_cell / 1000
    # end



    # procutil=Dict(Dagger.ThreadProc => 36.0)
    # @par cache=true procutil=procutil tseries =
    #     get_tseries(tseries_dir, tyh5_path, tseries_dset)

    @assign dag = (
        tseries_dir, fish_dir, zseries_dir, expName, recording_folder,
        fish_name, plot_dir, local_slm_dir, tyh5_path, test, tseries, zseries,
        zseries_xml, tseries_xml, zseries_zaxes, tseries_zaxis, result,
        exp_date, frame_rate, etl_vals, tseriesH, tseriesW, tseriesZ, tseriesT,
        zseriesH, zseriesW, zseriesZ, res, trial_order, slm_exp_dir,
        first_target_group, slm_num, power_per_cell, res2, slm1_power,
        slm2_power, voltageFile, res3, stimStartIdx, stimEndIdx, frameStartIdx,
        res4, target_groups, group_stim_freq, nStimuli, nTrials, zOffset,
        suite2p_dir, nwb_path, ndict, cell_traces, cell_masks, iscell, nCells,
        used_slm_dir, vol_rate, imaging2zseries_plane, zseries_xml
    )

    dag
end

function update_dag(r)
    @pun (
        tseries_dir, fish_dir, zseries_dir, expName, recording_folder,
        fish_name, plot_dir, slm_dir, tyh5_path, test, tseries, zseries,
        zseries_xml, tseries_xml, zseries_zaxes, tseries_zaxis, result,
        exp_date, frame_rate, etl_vals, tseriesH, tseriesW, tseriesZ, tseriesT,
        zseriesH, zseriesW, zseriesZ, res, trial_order, slm_exp_dir,
        first_target_group, slm_num, power_per_cell, res2, slm1_power,
        slm2_power, voltageFile, res3, stimStartIdx, stimEndIdx, frameStartIdx,
        res4, target_groups, group_stim_freq, nStimuli, nTrials, zOffset,
        suite2p_dir, nwb_path, ndict, cell_traces, cell_masks, iscell, nCells
    ) = r.thunks

    if ~in(:somekey, keys(r.thunks))
        somekey = @thunk get_somekey()
        @assign r.thunks = somekey
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

function Base.getindex(r::Recording, k)
    # collect(r.thunks[k])
    reify(r.thunks[k])
end

#### main analysis functions

"I thunk nothing."
ithunk() = thunk(identity)(nothing)

function is_setting(settings, key, continuation_thunk=ithunk())
    if key in keys(settings)
        settings[key]
    else
        reify(continuation_thunk)
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

function get_tseries(tseries_dir, tyh5_path, tseries_dset)
    if ~isnothing(tseries_dset)
        if isnothing(tyh5_path)
            tyh5_path = tseries_dir*".ty.h5"
        end
        h5, tseries = lazy_read_tyh5(tyh5_path, tseries_dset);
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