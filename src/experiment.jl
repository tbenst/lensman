using Lensman.DictThreadSafe
# import Lensman.DictThreadSafe: dictsrv
struct Recording
    uri::String
    settings::Dict # not thread safe
    # cache::DictSrv # thread safe
    cache::Dict # thread safe
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
function Recording(uri, settings)
    Recording(
        uri, merge(DEFAULT_SETTINGS, settings), Dict()
    )
end

function Base.getindex(r::Recording, k)
    if k in keys(r.cache)
        r.cache[k]
    else
        result = _compute_rec(r, k)
        r.cache[k] = result
        result
    end
end


Base.setindex(e::Recording, k, v) = (e.cache[k] = v)

"""Compute and assign local variables.
```
@experiment (etlVals, Z, H, W) = "2021-08-03/fish1/TSeries-1"
```

will expand to:
```
etlVals, Z, H, W = (
    compute("2021-08-03/fish1/TSeries-1", Val(:etlVals)),
    compute("2021-08-03/fish1/TSeries-1", Val(:Z)),
    compute("2021-08-03/fish1/TSeries-1", Val(:H)),
    compute("2021-08-03/fish1/TSeries-1", Val(:W)),
```
"""
macro experiment(ex)
    @assert ex.head == :(=) "no `=` found in expression."
    vars = ex.args[1] # :((etlVals, Z, H, W))
    exp_name = ex.args[2] # "2021-08-03/fish1/TSeries-1"
    vars_str = map(string, vars.args) # ["etlVals", "Z", "H", "W"]
    esc(:($vars = compute_rec($exp_name, $vars_str)))
    # ex
end

function compute_rec(r::Recording, ::Val{:tseries_dir})
    for root_dir in r.settings[:tseries_root_dirs]
        path = joinpath(root_dir, r.uri)
        if isdir(path)
            return path
        end
    end
    "could not find experiment"
end

function compute_rec(r::Recording, ::Val{:fish_dir})
    @pun tseries_dir = r
    sp = splitpath(tseries_dir)[1:end-1]
    joinpath(sp...)
end

function compute_rec(r::Recording, ::Val{:plot_dir})
    @pun fish_dir = r
    plot_dir = joinpath(fish_dir, r.settings[:rel_plot_dir])
end

function compute_rec(r::Recording, ::Val{:exp_name})
    @pun tseries_dir = r
    expName = splitpath(tseries_dir)[end]
end

function compute_rec(r::Recording, ::Val{:recording_folder})
    @pun tseries_dir = r
    recording_folder = splitpath(tseriesDir)[end-2]
end

function compute_rec(r::Recording, ::Val{:fish_name})
    @pun tseries_dir = r
    fish_name = splitpath(tseriesDir)[end-1]
end

function compute_rec(r::Recording, ::Val{:slm_dir})
    @pun fish_dir = r
    slm_dir = joinpath(fish_dir, "slm")
end

function compute_rec(r::Recording, ::Val{:tyh5_path})
    @pun tyh5_path = r.settings
    if isnothing(tyh5_path)
        @pun tseries_dir = r
        tyh5_path = tseries_dir*".ty.h5"
    end
    tyh5_path
end

function compute_rec(r::Recording, ::Val{:tseries})
    @pun (tyh5_path, tseries_dset) = r.settings
    @pun tseries_dir = r
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
