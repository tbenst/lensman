using Lensman.DictThreadSafe
using Dagger
import Dagger: @par
# import Lensman.DictThreadSafe: dictsrv
struct Recording
    uri::String
    thunks::Dict
    Recording(uri, settings) = new(
        uri, make_dag(uri, merge(DEFAULT_SETTINGS, settings)))
end

# Recording(uri, settings) = Recording(uri, settings, dictsrv(Dict()))
DEFAULT_SETTINGS = Dict(
    :tseries_root_dirs => [
        "/oak/stanford/groups/deissero/users/tyler/b115",
        "/scratch/b115",
        "/data/dlab/b115",
        "/mnt/deissero/users/tyler/b115"
    ],
    :rel_plot_dir => "plots",
    :tseries_dset => nothing,
    :tyh5_path => nothing,
)

"Master DAG for all computations."
function make_dag(uri, settings)
    dag = Dict()
    @pun (rel_plot_dir, tseries_dset) = settings
    @par cache=true begin
        tseries_dir = get_tseries_dir(uri, settings[:tseries_root_dirs])
        fish_dir = get_fish_dir(tseries_dir)
        expName = (x->splitpath(x)[end])(tseries_dir)
        recording_folder = (x->splitpath(x)[end-2])(tseries_dir)
        fish_name = (x->splitpath(x)[end-1])(tseries_dir)
        plot_dir = joinpath(fish_dir, rel_plot_dir)
        slm_dir = joinpath(fish_dir, "slm")
        tyh5_path = get_tyh5_path(settings[:tyh5_path], tseries_dir)

        tseries = get_tseries(tseries_dir, tyh5_path, tseries_dset)
    end
    @assign dag = (tseries_dir, fish_dir, tyh5_path, tseries)
    dag
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
    collect(r.thunks[k])
end

#### main analysis functions


function get_tseries_dir(uri, tseries_root_dirs)
    for root_dir in tseries_root_dirs
        path = joinpath(root_dir, uri)
        if isdir(path)
            return path
        end
    end
    @error "could not find experiment"
end

function get_fish_dir(tseries_dir)
    sp = splitpath(tseries_dir)[1:end-1]
    joinpath(sp...)
end

function get_tyh5_path(settings_tyh5_path, tseries_dir)
    if isnothing(settings_tyh5_path)
        tyh5_path = tseries_dir*".ty.h5"
    end
    tyh5_path
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

"Return a thunk or a value."
function _compute_rec(r::Recording, s::Symbol)
    result = compute_rec(r, Val(s))
    r.cache[s] = result
    result
end
