ENV["DISPLAY"] = "localhost:11"
##
using ImageView, Cairo
using Lensman, Images, Glob, StatsBase, ImageMagick, ImageTransformations,
    MicroscopyLabels, AxisArrays, ImageAxes, EndpointRanges,
    CoordinateTransformations, Rotations, FileIO
import Lensman: @pun, @assign
L = Lensman
##
resources = Resources();
uri = "2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061"
wt = Recordings[uri](;
    resources..., tseries_read_strategy=:hwzt,
    tyh5_path = "/data/dlab/b115/2021-06-01_wt-chrmine_h2b6s/fish4/TSeries-lrhab-control-118trial-061_red_kalman.h5",
);

## tif approach
wt_path = "2021-03-08_wt-chrmine-gc6f/fish1/SingleImage-1100nm-1024-037/" # may be unfair due to GC6f!
wtp2 = "2021-03-16_wt-chrmine_h2b6s/fish2/SingleImage-1115nm-512-027-new"
wtp3 = "2021-03-16_wt-chrmine_h2b6s/fish2/SingleImage-115nm-512-028"
wtp4 = "2021-02-16_h2b6s_wt-chrmine/fish3/SingleImage-1000nm-didnt-move-up-1024-037"
wtp5 = "2021-03-08_wt-chrmine-gc6f/fish1/SingleImage-1200nm-1024-031"
wtp6 = "2021-03-08_wt-chrmine-gc6f/fish1/SingleImage-1125nm-1024-045"
rs_path = "2021-05-18_rsChRmine_h2b6s_6dpf/fish3/SingleImage-1024-1115nm_habplane-061"
rs_path2 = "2021-03-16_rschrmine_h2b6s/fish3/SingleImage-115nm-512-028"
rs_path3 = "2021-05-18_rsChRmine_h2b6s_6dpf/fish5/SingleImage-1024-1115nm_30plane-063"
rs_path4 = "2021-05-18_rsChRmine_h2b6s_6dpf/fish1/SingleImage-1024-1115nm_habplane-055"
rsp5 = "2021-03-16_rschrmine_h2b6s/fish3/SingleImage-1115nm-512-027-new"

function load_img(rel_path, root_dirs, gamma=0.5)
    wt_tif_dir = L.find_folder(rel_path, root_dirs)
    @show wt_tif_dir
    fp = L.glob_one_file("*Ch2*.tif", wt_tif_dir)
    wt_im = ImageMagick.load(fp);
    @show mean(wt_im)
    wt_im = adjust_histogram(wt_im, GammaCorrection(gamma));
    wt_im = RGB.(imadjustintensity(wt_im));
    channelview(wt_im)[[2,3],:,:] .= 0
    wt_im
end

# wt_im = load_img(wt_path, wt[:tseries_root_dirs]); # 0.0012
# wt_im = load_img(wtp2, wt[:tseries_root_dirs]); # 1e-3, bad plane
# wt_im = load_img(wtp3, wt[:tseries_root_dirs]) # 9e-5, terrible
# wt_im = load_img(wtp4, wt[:tseries_root_dirs]) # 8e-4, forebrain plane / crooked
# wt_im = load_img(wtp5, wt[:tseries_root_dirs]) # 8e-4,
wt_im = load_img(wt_path, wt[:tseries_root_dirs],0.4) # 1e-3; looks great; gcamp contamination..?,
rs_im = load_img(rs_path, wt[:tseries_root_dirs],0.3); # 0.00042
# rs_im = load_img(rs_path2, wt[:tseries_root_dirs]); # 0.0013, bad plane
# rs_im = load_img(rs_path3, wt[:tseries_root_dirs]); # 5e-5, bad plane
# rs_im = load_img(rs_path4, wt[:tseries_root_dirs]); # 1e-3, grainy
# rs_im = load_img(rsp5, wt[:tseries_root_dirs]); # 1e-3, bad plane

W = size(rs_im,2)
lateral_unit = microscope_lateral_unit(W)
##
trfm = recenter(RotMatrix(-5*pi/256), center(rs_im));
rs_rot = warp(rs_im, trfm)
trfm = recenter(RotMatrix(-2*pi/256), center(wt_im));
wt_rot = warp(wt_im, trfm)
rs_img = AxisArray(rs_rot, (:y, :x), (lateral_unit,lateral_unit))
rs_img = rs_img[190:end-200,20:end-230]
wt_img = AxisArray(wt_rot, (:y, :x), (lateral_unit,lateral_unit))
wt_img = wt_img[180:end-190,20:end-220]
scalebar!(rs_img, 50u"μm")
scalebar!(wt_img, 50u"μm")
rs_wholeH = size(rs_img,2)
wt_wholeH = size(wt_img,2)
# imshow(rs_img)
# imshow(wt_img)
# mean(rs_img), mean(wt_img)
Array(wt_img)
Array(rs_img)
plot_dir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/1a"
pwt = joinpath(plot_dir, "wt-chrmine_whole-brain.png")
prs = joinpath(plot_dir, "rs-chrmine_whole-brain.png")
FileIO.save(pwt, wt_img)
FileIO.save(prs, rs_img)
println(pwt)
println(prs)
## zoomed in image
wt_im = load_img(wt_path, wt[:tseries_root_dirs],1) # 1e-3; looks great; gcamp contamination..?,
channelview(wt_im)[[2,3],:,:] .= 0
rs_im = load_img(rs_path, wt[:tseries_root_dirs],1); # 0.00042
channelview(rs_im)[[2,3],:,:] .= 0
rs_im
trfm = recenter(RotMatrix(-5*pi/256), center(rs_im));
rs_rot = warp(rs_im, trfm)
trfm = recenter(RotMatrix(-2*pi/256), center(wt_im));
wt_rot = warp(wt_im, trfm)
rs_img = AxisArray(rs_rot, (:y, :x), (lateral_unit,lateral_unit))
rs_img = rs_img[190:end-200,20:end-230]
wt_img = AxisArray(wt_rot, (:y, :x), (lateral_unit,lateral_unit))
wt_img = wt_img[180:end-190,20:end-220]
rs_img = adjust_histogram(rs_img, GammaCorrection(0.5));
rs_img = adjust_histogram(rs_img, ContrastStretching(t=0.2,slope=1.5))
channelview(rs_img)[[2,3],:,:] .= 0
cmax = quantile(channelview(rs_img[:])[1,:],0.999)
rs_img = Array(imadjustintensity(rs_img,(0,cmax)))
wt_img = adjust_histogram(wt_img, GammaCorrection(1.0));
# wt_img = adjust_histogram(wt_img, ContrastStretching(t=0.2,slope=1.5))
channelview(wt_img)[[2,3],:,:] .= 0
cmax = quantile(channelview(wt_img[:])[1,:],0.999)
wt_img = Array(imadjustintensity(wt_img,(0,cmax)))
# imshow(rs_img)
# imshow(wt_img)
n = 70
# h = 388
h = 360
# w = 587
w = 465
rs_img = rs_img[h:h+n,w:w+n]
channelview(rs_img)[[2,3],:,:] .= 0
curH = size(rs_img,1)
uscale = rs_wholeH/curH
rs_img = imresize(rs_img, W,W)
rs_img = AxisArray(rs_img, (:y, :x), (lateral_unit,lateral_unit) ./ uscale)
scalebar!(rs_img, 5u"μm")
Array(rs_img)

n = 70
# h = 388
h = 230
# w = 587
w = 410
wt_img = wt_img[h:h+n,w:w+n]
channelview(wt_img)[[2,3],:,:] .= 0
@show curH = size(wt_img,1)
uscale = wt_wholeH/curH
wt_img = imresize(wt_img, W,W)
wt_img = AxisArray(wt_img, (:y, :x), (lateral_unit,lateral_unit) ./ uscale)
scalebar!(wt_img, 5u"μm")

pwt = joinpath(plot_dir, "wt-chrmine_closeup.png")
prs = joinpath(plot_dir, "rs-chrmine_closeup.png")
FileIO.save(pwt, wt_img)
FileIO.save(prs, rs_img)
println(pwt)
println(prs)
