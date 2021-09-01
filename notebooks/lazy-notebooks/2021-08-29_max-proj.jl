ENV["DISPLAY"] = "localhost:11"
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
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    ](
    ;resource...,
    # tseries_dset=nothing, lazy_tiff=false
);
##
@pun (tseries, tseriesT, fish_dir, exp_name) = r;
##
size(tseries)
# size(tseries.dset)
##
tseries_subset = tseries[:,:,:,1:20:tseriesT];
##
maxTproj = mapslices(x->quantile(x[:],0.90), tseries_subset,dims=4);
shiftedmaxTproj = shift_every_other_row(maxTproj,-3);
imshow(shiftedmaxTproj)

##
im = shiftedmaxTproj;
im = Gray.(imadjustintensity(im));
# im = Gray.(imadjustintensity(im, (otsu_threshold(im),maximum(im))))
# im = convert(Array{Gray{N0f16}}, im);
im = im[:,:,:,1];
Gray.(im[:,:,5])
path = joinpath(fish_dir, exp_name*"lstm_stack_float64.tif")
save(joinpath(fish_dir, exp_name*"lstm_stack_float64.tif"), im);
path
# imshow(im)
# ImageMagick.save("/tmp/lstm_stack_uint16_gamma0.5.tiff", im);