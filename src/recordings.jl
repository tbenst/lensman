# TODO: we should have a Fish struct so we can be more DRY?
# also, are we just reinventing OOP here...? but a crummy version where hard to extend?
# => rebuttal: at least our methods are more modular, and pure, and better composition

"return a callable that can add nodes to DAG with kwargs"
function modifiable_recording(uri ;nodes...)
    uri => (;new_nodes...) -> update_recording_dag(DAG(;
        merge(DEFAULT_RECORDING_SETTINGS, nodes, new_nodes,
        Dict(:uri => uri))...))
end

DEFAULT_RECORDING_SETTINGS = Dict(
    :tseries_root_dirs => [
        "/oak/stanford/groups/deissero/users/tyler/b115",
        "/scratch/b115",
        "/scratch/b118",
        "/scratch2",
        "/data/dlab/b115",
        "/data/dlab/b118",
        "/data/b115",
        "/mnt/deissero/users/tyler/b115",
        "/mnt/deissero/users/tyler/b118"
    ],
    :slm_root_dirs => [
        "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files",
        "/mnt/deissero/users/tyler/b115/SLM_files",
        "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/",
        "/mnt/deissero/users/tyler/b115/SLM_files/",
        "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/",
        "/data/b115/SLM_files",
        "/data/dlab/b118/SLM_files"
    ],
    :rel_analysis_dir => "raw",
    :tseries_dset => nothing, # if nothing, use raw tiffs
    :tyh5_path => nothing,
    :slm_dir => "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/",
    # TODO: should we use `missing` instead of these strings..?
    :oir_dir => "OIR_DIR_NOT_DEFINED",
    :suite2p_dir => nothing,
    :warp_prefix => "WARP_PREFIX_NOT_DEFINED",
    # :tyh5, :lazy_tyh5, :lazy_hwzt, :lazy_tiff, or :tiff ("else")
    :tseries_read_strategy => :tiff,
    :region_mask_strategy => :cmtk,
    :lazy_tiff => true,
    :window_secs => 5,
    :cells_df_f_win_secs => 2,
    :cells_df_f_padding => 2,
    :ants_no_run => true, # TODO: make false once tested
    :oir_920_name => "OIR_920_NAME_NOT_DEFINED",
    :oir_820_name => "OIR_820_NAME_NOT_DEFINED",
    :zbrain_warp_prefix => "ZBRAIN_WARP_PREFIX_NOT_DEFINED",
    :mm_warp_prefix => "",
    :rostral => :right,
    :dorsal => :up,
    :notes => "",
    :zseries_name => missing,
    :genotype => "", # not yet used
    :registration_type => :dont_run,
    :roi_method => :suite2p, # or :cellpose
    :channel_str => "Ch3",
    :tseries_suffix => "raw", # or "lstm" or "kalman"
    :trial_average_path => nothing,
    :trial_average_dset => nothing
    # slm_dir => "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files",
    # slm_dir => "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files",
    # slm_dir => "/mnt/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment",
)

struct RecordingsWrapper
    dict::Dict
    RecordingsWrapper(args...;kwargs...) = new(Dict(args...;kwargs...))
end

Base.getindex(r::RecordingsWrapper, uri) = uri in keys(r.dict) ? r.dict[uri] : modifiable_recording(uri)[2]

# set tseries_dset = nothing to force use of tiff files
Recordings = RecordingsWrapper(
    modifiable_recording(
        "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061";
        zseries_name="ZSeries-structural-840nm-057"
    ),
    modifiable_recording(
        "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061";
        tseries_dset = "/imaging/PerVoxelLST1M_actually_shared-separate_bias_hidden_init_from_pretrained-2021-06-21_6pm",
        tyh5_path="/scratch/b115/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061.ty.h5",
        tseries_read_strategy=:lazy_tyh5
    ),
    modifiable_recording(
        "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-118trial-122";
        zseries_name="ZSeries-structural-840nm-058",
        tseries_dset = "/imaging/LSTM_per-voxel-state_divide8192-2021-07-02",
        tyh5_path="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-118trial-122.ty.h5",
        suite2p_dir="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/s2p/suite2p",
        tseries_read_strategy=:lazy_tyh5,
        notes="""
            Only did one stim pulse??
        """
    ),
    modifiable_recording(
        "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123";
        zseries_name="ZSeries-structural-840nm-058",
        tseries_dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm",
        tyh5_path="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-1232021-06-21_6pm.ty.h5",
        suite2p_dir="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/s2p/suite2p",
        tseries_read_strategy=:lazy_tyh5
    ),
    modifiable_recording(
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062";
        zseries_name="ZSeries-structural-840nm-048",
        # tseries_dset = "/imaging/LSTM_per-voxel-state_divide8192-2021-07-02",
        tyh5_path="/scratch/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062lstm_divide8192_avgStim.h5",
        oir_dir="/data/dlab/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fishfrom_2021-07-13_rschrmine_h2b6s/fish1",
        oir_920_name="multimap_zseries_920nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom.oir",
        oir_820_name="multimap_zseries_820nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom-take2.oir",
        tseries_read_strategy=:lazy_tyh5,
        mm_warp_prefix="20210829T042045359"
        # mm_warp_prefix="20210729T154219098Z"
    ),
    modifiable_recording(
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061";
        zseries_name="ZSeries-structural-840nm-048",
        # tseries_dset = "/imaging/LSTM_per-voxel-state_divide2048-2021-07-02",
        # tseries_read_strategy=:lazy_tyh5,
        lazy_tiff=true,
        tyh5_path="/data/dlab/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061.ty.h5",
        oir_dir="/data/dlab/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fishfrom_2021-07-13_rschrmine_h2b6s/fish1",
        oir_920_name="multimap_zseries_920nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom.oir",
        oir_820_name="multimap_zseries_820nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom-take2.oir",
        mm_warp_prefix="20210829T042045359"
        # mm_warp_prefix="20210729T154219098Z"
        # zbrain_warp_prefix="20210727T204656077Z" # TODO: is this zbrain or for multimap...?
    ),
    modifiable_recording(
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069";
        zseries_name="ZSeries-structural-840nm-056",
        lazy_tiff=true,
        tyh5_path="",
        oir_dir="/data/dlab/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fishfrom_2021-07-13_rschrmine_h2b6s/fish2",
        oir_920_name="multimap_zseries_920nm_ch4-gad405_ch1-sert647_chr-gcamp-2xzoom-day2.oir",
        oir_820_name="multimap_zseries_820nm_ch4-gad405_ch1-sert647_chr-gcamp-2xzoom-day2_0001.oir",
        mm_warp_prefix="20210828T215549354",
        notes="""Originally, we warped on the 920 stack, but there was drift between that and
        820 (due to water under dish, argh..). so instead we warp on the 820 stack
        (prefix: 20210828T215549354)
        """
    ),
    modifiable_recording(
        "2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-046";
        tseries_dset="/imaging/raw"
    ),
    modifiable_recording(
        "2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-5ppc-048";
        tseries_dset="/imaging/raw"
    ),
    modifiable_recording(
        "2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-047";
        tseries_dset="/imaging/raw"
    ),
    modifiable_recording(
        "2021-01-19_chrmine_kv2.1_h2b6s_7dpf/fish2_nochrmine/TSeries-1024cell-32concurrent-4power-048";
    tseries_dset="/imaging/raw"
    ),
    modifiable_recording(
        "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-titration-192trial-061";
        zseries_name="ZSeries-structural-056"
    ),
    modifiable_recording(
        "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072";
        tyh5_path="/data/dlab/b115/2021-06-01_rsChRmine_h2b6s/fish3/TSeries-IPNraphe-118trial-072_kalman.h5",
        dset_name="kalman",
        tseries_read_strategy = :lazy_hwzt,
        zseries_name="ZSeries-structural-056",
        notes="actually from 2021-06-02!!"
    ),
    modifiable_recording(
        "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072";
        tyh5_path="/data/dlab/b115/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072_kalman.h5",
        dset_name="kalman",
        tseries_read_strategy = :lazy_hwzt,
        zseries_name="ZSeries-structural-056",
        notes="copied from 2021-06-01/fish3, where we accidentally saved."
    ),
    modifiable_recording(
        "2022-03-11_rschrmine_h2b6s_9dpf/fish1/TSeries-raphe-thalamus-118trial-002";
        zseries_name="ZSeries-structural-056"
    ),
    modifiable_recording(
        "2022-03-31_rschrmine_h2b6s_6dpf/fish1/TSeries-raphe-thalamus-118trial_no-autocorrect-350p-005";
        zseries_name="ZSeries-860nm_galvo-012"
    ),
    modifiable_recording(
        "2022-05-19_rschrmine_h2b6s_6dpf/fish1/TSeries-raphe-thalamus-5power-226trial-003";
        tseries_read_strategy=:lazy_tiff,
        zseries_name="ZSeries-structural-860nm-resonant-009"
    ),

)