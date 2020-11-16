##
true
##
# ENV["DISPLAY"] = "localhost:11.0" # remote-ssh
ENV["DISPLAY"] = ":0"
using FileIO, RegisterQD, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, Statistics, Random,
    ImageDraw, Lensman, ImageFiltering
import Base.Threads.@threads
using Unitful: μm, s
sio = pyimport("scipy.io")

##
# dataDir = "/mnt/deissero/users/tyler/b115/2020-10-12_elavl3-chrmine-kv2.1_h2b-6s_11dpf/zstack-1x/SingleImage-820_128avg-032"
dataDir = "/mnt/deissero/users/tyler/b115/2020-10-12_elavl3-chrmine-kv2.1_h2b-6s_11dpf/zstack-1.5x/ZSeries-09172019-1403-029"
# tif = joinpath(dataDir, "SingleImage-820_128avg-032_Cycle00001_Ch3_000001.ome.tif")
tif = joinpath(dataDir, "ZSeries-09172019-1403-029_Cycle00001_Ch3_000050.ome.tif")
##
img = ImageMagick.load(tif)
eq = ImageContrastAdjustment.AdaptiveEqualization(nbins = 64, rblocks = 8, cblocks = 8, clip = 0.1)
img = imadjustintensity(img)
# img = adjust_histogram(img, GammaCorrection(gamma = 0.5))
# eq = Equalization(nbins=64)
img = adjust_histogram(img, eq)
# img = adjust_histogram(adjust_histogram(img, GammaCorrection(gamma = 0.5)), eq)
# imadjustintensity(adjust_histogram(img, GammaCorrection(gamma = 0.4)))
# adjust_histogram(Gray.(img), Equalization())
# imadjustintensity(img)
# img = mapwindow(median,imadjustintensity(img), (3,3));
img
##


## READ UNITS AND MAKE NEURON TEMPLATE
microscope_units = read_microns_per_pixel(datadir)
## CREATE NUCLEUS TEMPLATE
nucleus_r = 2.75μm/microscope_units[1]
space = 0.3μm # WARN: need to adjust CirclePointRadius, too
cell_r = nucleus_r + space
template_size = Int16(ceil(2*(cell_r)))
# ensure odd size
template_size += abs(template_size % 2 - 1)
# we make large and later resize for anti-aliasing
nucleus_template = Gray.(zeros(513,513))
draw!(nucleus_template, Ellipse(CirclePointRadius(256,256,192)))
sz = (template_size, template_size)
σ = map((o,n)->0.75*o/n, size(nucleus_template), sz)
kern = KernelFactors.gaussian(σ)

nucleus_template = imresize(imfilter(nucleus_template, kern, NA()), sz)
##
# rescaled_img = adjust_gamma(imadjustintensity(img), 0.45)
rescaled_img = imadjustintensity(img)


# outsu mask
threshold = convert(Float64,otsu_threshold(rescaled_img))
threshold *= 1.5
fg_mask = rescaled_img.>threshold

image_fg = rescaled_img .* fg_mask

# imshow(image_fg.>0)


filt = imfilter(Float64.(image_fg), Float64.(nucleus_template), Fill(0.))
filt = Gray.((filt.-minimum(filt))./(maximum(filt).-minimum(filt)))

seeds = findlocalmaxima(filt)
seed_mask = zeros(Bool, size(filt))
seed_mask[seeds] .= true
seed_mask = label_components(seed_mask)
segments = watershed(Float64.(image_fg), seed_mask)


imshow(map(i->get_random_color(i), labels_map(segments)) .* (1 .-fg_mask)) 

##
seed_mask = zeros(Bool, size(filt))
seed_mask[seeds] .= true
imshow(Gray.(seed_mask))
