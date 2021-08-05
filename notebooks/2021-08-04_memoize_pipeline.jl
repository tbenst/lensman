using Lensman
import Lensman: @pun
using Dagger
import Dagger: @par
##
settings = Dict(
    :tseries_root_dirs => [
        "/data/b115"
    ],
    # :rel_plot_dir => "plots_denoised",
    # :tseries_dset => "/imaging/LSTM_per-voxel-state_divide2048-2021-07-02",
)

r = Recording(
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123",
    settings
)
@pun (tseries, tseries_dir, fish_dir, tyh5_path) = r
