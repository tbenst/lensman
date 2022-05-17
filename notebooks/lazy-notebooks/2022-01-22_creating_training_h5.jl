##
using Glob, NPZ, PyCall, Images
using Lensman, Thunks
import Lensman: @pun
L = Lensman
np = pyimport("numpy")
# import PyPlot
# plt = PyPlot
# matplotlib = plt.matplotlib

##
suite2p_dir = "/mnt/stanford_oak_projects/babelfish/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/s2p/suite2p"

## Test out some reading
plane_dirs = glob("plane*", suite2p_dir);
n_planes = length(plane_dirs);

iscell = npzread(joinpath(suite2p_dir, "combined", "iscell.npy"));
n_cells = size(iscell, 1);

##
# Create combined_stat by combining across planes (seems this is still necessary according to:
# https://github.com/MouseLand/suite2p/issues/655)
combined_stat = PyObject[]
# Assuming that plane files are indexed 0, 1, ..., n_planes-1
for z in 1:n_planes
    stat = np.load(joinpath(suite2p_dir, "plane$(z - 1)", "stat.npy"), allow_pickle = true)
    for k in keys(stat)
        # Use 1 indexing (Julia) rather than 0 indexing (suite2p/python)
        stat[k].update(iplane = z)  # can't do [] assignment in PyCall
    end
    combined_stat = vcat(combined_stat, stat)
end

# Pixel order for center is (y, x)
cell_centers = map(x -> np.array(get(x, "med")), combined_stat);
planes = map(x -> get(x, "iplane"), combined_stat);
for (cell_center, plane_index) in zip(cell_centers, planes)
    append!(cell_center, plane_index)
end

##
# Read sample trace and see what it looks like!
# z = 2
# F = np.load(joinpath(suite2p_dir, "plane$(z - 1)", "F.npy"), allow_pickle = true);
F_combined = np.load(joinpath(suite2p_dir, "combined", "F.npy"), allow_pickle = true);

using Plots
plot(F_combined[1, :])

##
# So, there was this region_masks.h5. Seems this is what we want (just need to make sure to
# remove None's, etc.)
# And then this pyarrow file (this one for kalman, and for this sequence specifically)
# using Arrow
# table_file = joinpath(mask_dir, "TSeries-lrhab-118trial-061_kalman.h5/kalman", "regions_df.arrow")
# table = Arrow.Table(table_file)
using HDF5


mask_dir = "/scratch/allan/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/h5_output"
println("Will read masks from: '$mask_dir'")
region_masks_h5_path = joinpath(mask_dir, "warped_region_masks.h5")
region_masks_h5 = h5open(region_masks_h5_path, "r", swmr = true)

##
# Loaded data to see which specific tseries it is. Seems it's this one, length 7056
resources = Resources(zbrain_dir = "/data/b115/atlas/zbrain");
@pun (zbrain_dir, zbrain_masks) = resources;

recording_name = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
options = Dict(
    :rel_analysis_dir => "plots", :lazy_tyh5 => true, :lazy_tiff => true,  # changed second
    :tseries_root_dirs => ["/data/b115"], :slm_root_dirs => ["/data/b115/SLM_files"]
)
recording = Recordings[recording_name](; options..., resources...);
@pun (imaging2zseries_plane, window_len,
    tseries, tseriesH, tseriesW, tseriesZ, trial_order,
    zbrain_mask_names, target_groups,
    stim_start_idx, stim_end_idx) = recording;

##
# Try to remove things that are outside of brain (not sure if this means just
# not in any region)
region_names = zbrain_masks["MaskDatabaseNames"][1:end-2]

##
# Read all masks by name, could change `zbrain_mask_names` from Recording to `region_names` from Resource
masks = map(n -> L.try_read_mask(region_masks_h5, zbrain_mask_names, imaging2zseries_plane, n),
    region_names);

##
# Remove masks returning as `nothing` (assuming correct flushing when creating region_masks.h5, this should
# just be masks saved as 0 because no pixels overlapping with imaged - should only be a few)
valid_mask_indices = (masks .!== nothing)
println("Dropping $(length(masks) - sum(valid_mask_indices)) valid masks out of $(length(masks)), because = nothing")
region_names = region_names[valid_mask_indices];
masks = masks[valid_mask_indices];

##
# Remove masks with no voxels imaged
valid_mask_indices = map(m -> sum(m) > 0, masks)
println("Dropping a further $(length(masks) - sum(valid_mask_indices)) valid masks out of $(length(masks)), because 0 in imaged plane")
region_names = region_names[valid_mask_indices];
masks = masks[valid_mask_indices];

##
# Left vs. right can be added, just copy code over

##
masks = map(m -> imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);

##
# Visualize sum of all masks to see what it looks like
brain_mask = sum(masks)
img1 = brain_mask[:, :, 5]
Gray.(convert(Matrix{Float64}, img1))

##
# See how many neuron centers should be discarded

region_name_to_neurons = Dict(k => Vector{Int}() for k in region_names)

cell_in_brain = ones(Bool, length(cell_centers))
println("Number of cells at start: $(sum(cell_in_brain))")

for (cell_number, cell_center) in enumerate(cell_centers)
    cc = round.(Int, cell_center)
    if brain_mask[cc[1], cc[2], cc[3]] == 0
        cell_in_brain[cell_number] = false
    else
        for (region_name, mask) in zip(region_names, masks)
            # if mask[Int(round(y)), Int(round(x)), plane_index]
            if mask[cc[1], cc[2], cc[3]]
                append!(region_name_to_neurons[region_name], cell_number)
            end
        end
    end
end

println("Found $(sum(cell_in_brain)) out of $(length(cell_in_brain)) cells to be in brain")

##
cell_centers = Matrix{Float64}(reduce(hcat, cell_centers)')  # not sure why this cast is needed

##
# Can drop file size by converting to Float32 or Int where possible

##
println("cell_centers has shape $(size(cell_centers))")
h5open("example2.h5", "w") do fid
    fid["fluorescence"] = F_combined
    println(typeof(F_combined))
    fid["cell_in_brain"] = cell_in_brain
    println(typeof(cell_in_brain))
    fid["iscell"] = iscell
    println(typeof(iscell))
    fid["centroids"] = cell_centers  # (y, x, z), where z is just plane index
    println(typeof(cell_centers))


    g = create_group(fid, "regions")
    for (region_name, cell_indices) in region_name_to_neurons
        g[region_name] = cell_indices
    end
end

##
# Maybe try loading and making sure data is there
using HDF5
using Plots
h5open("example2.h5", "r") do fid
    x = read(fid["cell_in_brain"])
    println("cell_in_brain has shape: $(size(x)) and type $(typeof(x))")
    println("Number of cells in brain: $(sum(x))")
    x = read(fid["centroids"])
    println("centroids has shape: $(size(x)) and type $(typeof(x))")
    x = read(fid["fluorescence"])
    println("fluorescence has shape: $(size(x)) and type $(typeof(x))")
    plot(x[2, :])
    x = read(fid["iscell"])
    println("iscell has shape: $(size(x)) and type $(typeof(x))")
    println(x)
    # region_name = "Diencephalon -"
    # g = fid["regions"]
    # x = read(g[region_name])
    # println("Number of regions with cells in them: $(length(keys(g)))")
    # println("regions/$region_name has shape: $(size(x)) and type $(typeof(x))")
    # println("region '$region_name' has $(length(x)) cells in it")
end
# Proper way to do this is to construct mask for each neuron, and keep n_pixels/percent
# overlap with each region
# Also, these regions are not mutually exclusive so keeping all of them.