Resources(;nodes...) = update_resource_dag(
    DAG(;merge(DEFAULT_RESOURCE_SETTINGS, nodes)...))

# Recording(uri, settings) = Recording(uri, settings, dictsrv(Dict()))
DEFAULT_RESOURCE_SETTINGS = Dict(
    :zbrain_dir => "/mnt/deissero/users/tyler/zbrain",
    :zbrain_units => (0.798μm, 0.798μm, 2μm),
)

"Master DAG for all resources. These are shared across recordings."
function update_resource_dag(resource)
    @pun (zbrain_dir,) = resource
    @lazy begin
        h2b_zbrain = read_zbrain_line("$zbrain_dir/AnatomyLabelDatabase.hdf5",
            "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish")
    end

    @assign resource = h2b_zbrain
    resource
end
