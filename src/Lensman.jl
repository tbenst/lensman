module Lensman

using AxisArrays, ANTsRegistration, NIfTI, ImageMagick, Images,
    ImageDraw, ImageFiltering, PyCall, MAT, Dates, Unitful,
    # DataStructures,
    Statistics, SharedArrays, CSV, DataFrames, Suppressor, Plots,
    LinearAlgebra, LibExpat, LightXML, RollingFunctions, HypothesisTests,
    EllipsisNotation, HDF5, Distributed, Thunks, NPZ, SparseArrays,
    H5Sparse, Folds, Arrow, ProgressMeter, JuMP, ReusePatterns, StatsBase
import Base.Threads: @threads, @spawn, @sync
using Distributed
import Unitful:μm
import Base.Threads.@threads
sio = pyimport("scipy.io")
# using Dagger
# https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules

include("Types.jl")
include("generic.jl")
include("algorithms.jl")
include("files.jl")
include("plots.jl")
include("magic_numbers.jl")
include("hdf5_threads.jl")
include("dict_threads.jl")
include("images.jl")
include("distributed.jl")
include("dataframes.jl")
include("Bruker.jl")
include("core.jl")
include("tidying.jl")
include("resources.jl")
include("recording_dag.jl")
include("recordings.jl")
include("stats.jl")


export read_microns_per_pixel,
    read_mask,
    zbrain_units,
    antsApplyTransforms,
    apply_cmtk_transform,
    save_region_masks_cmtk,
    getMaskNameIncludes,
    read_zseries,
    create_nuclei_mask,
    findNeurons,
    create_slm_stim,
    mask2ind,
    readZseriesTiffDir,
    cartesianIndToArray,
    arrCart2arr,
    findNeuronsAaron,
    cartIdx2SeanTarget,
    cartIdx2Array,
    stripLeadingSpace,
    findTTLStarts,
    findTTLEnds,
    readTseriesTiffDir,
    getFramePlane,
    tseriesTiffDirMetadata,
    addTargetsToImage,
    getStimTimesFromVoltages,
    getTimeFromFilename,
    getStimFrameIdx,
    getExpData,
    getTrialOrder,
    loadTseries,
    write_experiment_to_tyh5,
    get_slm_stim_masks,
    init_workers,
    trialAverage,
    colorImage,
    countPulsesMaxGap,
    findMatGroups,
    allWithin,
    rectangle,
    extractTrace,
    imageJkalmanFilter,
    stonerStimGroups,
    rpadArray,
    calc_concurrency_score,
    randomSwaps,
    minMaxIntersection,
    count_concurrency,
    noDups,
    printGroups,
    aa2a,
    calcGroupsPerCell,
    stonerPerm,
    findIdxOfClosestElem,
    mapTargetGroupsToPlane,
    constructGroupMasks,
    imageCorrWithMask,
    Green,
    sym_adjust_gamma,
    zeroAdjust,
    perm,
    cartIdx2SeanTarget512,
    reconstructA,
    addCellIdx,
    f_lookup_cellidx,
    entangledInfluenceMaps,
    makeCellsDF,
    constructROImasks,
    plotStim,
    printMatKeys,
    getImagingPockels,
    parseXML,
    createMarkPointElement,
    medianfilt,
    read_gpl,
    max_95ci,
    min_95ci,
    getzoffset,
    spiral_size,
    microscope_lateral_unit,
    getMatStimFreq,
    getSLMnum,
    slmpower,
    # write_markpoints,
    write_trial_order,
    regex_glob,
    loadBOT,
    markpoints_magic_numbers,
    read_markpoints_series,
    write_markpoints,
    antsApplyTransforms,
    vecvec2mat,
    read_markpoint_groups,
    balanced_transition_order,
    lazy_read_tyh5,
    read_tyh5,
    read_oir_units,
    ants_register,
    avgTseries,
    match_histogram,
    mutual_information,
    glob_one_file,
    transparent_cmap,
    funprod,
    zeroToOne,
    read_xml,
    read_first_zaxis,
    read_all_zaxis,
    ez_gamma,
    timeseries_df_f,
    Recording, compute_rec, debug,
    shift_every_other_row, get_random_color,
    df_f_per_voxel_per_trial_from_h5,
    get_df_f_per_trial_dataframe,
    stim_roi_threshold,
    Recordings,
    trial_average_from_h5,
    Resources,
    read_registered_mask,
    influence_map,
    LazyTy5,
    noop, LazyTiff,
    opening_median, opening_nonbool
end