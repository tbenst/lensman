ENV["DISPLAY"] = "localhost:12"
##
using ImageView
using Lensman, Images, StatsBase
import Lensman: @pun, @assign
L = Lensman
##
resource = Resources();
r = Recordings[
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    # "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
    ](
    ;resource...,
    tseries_read_strategy = :lazy_hwzt,
    tyh5_path=nothing, tseries_dset=nothing
);
@pun (tseries_dir, tyh5_path) = r
##
@pun (tseries, tseriesT, fish_dir, exp_name) = r;
typeof(tseries)
##
size(tseries)
# size(tseries.dset)
##
tseries_subset = tseries[:,:,:,1:20:tseriesT];
##
maxTproj = mapslices(x->quantile(x[:],0.90), tseries_subset,dims=4);
# only for "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
# shiftedmaxTproj = shift_every_other_row(maxTproj,-3);
shiftedmaxTproj = maxTproj
imshow(shiftedmaxTproj)

##
im = shiftedmaxTproj;
im = Gray.(imadjustintensity(im));
##
# im = Gray.(imadjustintensity(im, (otsu_threshold(im),maximum(im))))
# im = convert(Array{Gray{N0f16}}, im);
im = im[:,:,:,1];
Gray.(im[:,:,5])
path = joinpath(fish_dir, exp_name*"kalman_stack_float64.tif")
save(joinpath(fish_dir, exp_name*"kalman_stack_float64.tif"), im);
path
##
maxZproj = maximum(im, dims=3)[:,:,1]
maxZproj = adjust_histogram(maxZproj, GammaCorrection(0.8))
im2 = maxZproj[175:310,100:235]
save(joinpath(fish_dir, exp_name*"closeup.png"),im2)
im2