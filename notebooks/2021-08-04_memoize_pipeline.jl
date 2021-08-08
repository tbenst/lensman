using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, Random
import Lensman: @pun, @assign
using Thunks
np = pyimport("numpy")
L = Lensman
# using Dagger
# import Dagger: @par
# import Dagger.Sch: ThunkOptions
##

r = Recording(
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123";
    zseries_name="ZSeries-structural-840nm-058",
    # tseries_dset = "/imaging/LSTM_per-voxel-state_divide2048-2021-07-02",
    tseries_dset = "/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm",
    tyh5_path="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-1232021-06-21_6pm.ty.h5",
    suite2p_dir="/scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/s2p/suite2p"
    # tseries_root_dirs = [
    #     "/data/b115"
    # ],
);
##
@pun (nCells, cell_masks, tseries_dir, suite2p_dir) = r
##
# @pun (tyh5_path, tseries_dir, fish_dir, test, exp_date, frame_rate, etl_vals) = r;
##
# @pun (stimStartIdx, stimEndIdx, tseries, tseriesH, tseriesW, tseriesZ) = r;
# size(tseries)
##
avg_tseries = avgTseries(tseries_dir);


##
labels = zeros(UInt16, size(avg_tseries)...)
# TODO: find non unique XYZ tuples & keep row with largest weight.

for c in unique(cell_masks.cell_id)
    pix_mat = Array(cell_masks[cell_masks.cell_id .== c,[:y,:x,:z]]) .+ 1
    voxel_mask = mapslices(x->CartesianIndex(x...),pix_mat,dims=2)
    labels[voxel_mask] .= c
end
##
colored_neurons = map(get_random_color, (labels[:,:,6]));
colored_neurons
shift_every_other_row(colored_neurons,-3)

##


## example for extending DAG....
# @thunk begin
#     suite2p_dir = L.is_setting(:suite2p_dir, get_suite2p_dir(r.thunks[:tseries_dir]))
# end
##
suite2p_dir = thunk(L.is_setting)(r.settings, :suite2p_dir, thunk(L.get_suite2p_dir)(r.thunks[:tseries_dir]))
@assign r.thunks = suite2p_dir
@pun suite2p_dir = r

## pattern for extending DAG (ugly! ergonomics should be improved..macro?)
if ~in(:somekey, keys(r.thunks))
    somekey = @thunk get_somekey()
    @assign r.thunks = somekey
end