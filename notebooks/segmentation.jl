##
ENV["DISPLAY"] = "localhost:11.0"
using ImageView, Images, ImageView, ImageSegmentation, ImageMagick, Random,
    ImageDraw, Lensman, ImageFiltering
import Base.Threads.@threads
using Unitful: μm, s
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end
##
datadir = "/data/dlab/zfish_bruker/2020-08-04_elavl3-ChRmine-Kv2.1_elavl3-H2B-GCaMP6s/fish1/fov1-hab-brainwide/green-003/"
imgpath = joinpath(datadir,"green-003_Cycle00001_Ch3_000001.ome.tif")
img = ImageMagick.load(imgpath);
##
segment_nuclei(img)
##
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
##
# save("/home/tyler/Desktop/zfish.png", adjust_gamma(rescaled_img,1.5))
##
rescaled_img = imadjustintensity(img)
threshold = convert(Float64,otsu_threshold(rescaled_img))
threshold *= 1.2
bw = rescaled_img.>threshold
Gray.(bw)
imshow(rescaled_img)
##
seeds = findlocalmaxima(filt)
seed_mask = zeros(Bool, size(filt))
seed_mask[seeds] .= true
seed_mask = label_components(seed_mask)
segments = watershed(Float64.(rescaled_img), seed_mask)
imshow(map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw)) 
##
dist = 1 .- distance_transform(feature_transform(bw));
markers = label_components(dist .< -15);
segments = watershed(dist, markers)
imshow(map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw))
# histeq(img,1024)
# histeq(adjust_gamma(img, 1),1024)
# img
##

##
Gray.(rescaled_img.>threshold)


##
using PyCall
py"""
import os
import numpy as np
import scipy
from scipy import ndimage
import skimage
from skimage.morphology import disk, watershed
from skimage.filters import rank
from skimage.filters.rank import otsu, autolevel_percentile, enhance_contrast_percentile
from skimage.filters import threshold_otsu
from skimage.feature import peak_local_max
from skimage.exposure import adjust_gamma
def segment_nuclei(img):
    maxval = np.max(img)
    img_gamma = np.power(img,0.45) * maxval / np.power(maxval,0.45)
    img_gamma = img_gamma.astype('uint16')

    # use image histogram to filter out background
    nbins = 200
    hist, bin_edges = np.histogram(np.ravel(img_gamma), bins=nbins)
    hthres = bin_edges[np.argmax(hist) + 6]
    #print "Histogram threshold: " + str(float(hthres) / np.max(img_gamma))
    fg_mask = (img_gamma >= hthres)

    # autolevel_percentile
    selem = disk(4)
    img_level = autolevel_percentile(img_gamma, selem=selem, p0=.1, p1=.90)
    img_level *= fg_mask  # apply otsu mask
    
    # otsu threshold 2 to binary image
    gthres = threshold_otsu(img_level)
    gthres*=1.2
    img_gthres = img_level >= gthres
    
    # find cell centers
    kern = disk(2)
    distance = ndimage.distance_transform_edt(img_gthres)
    feat_img = skimage.feature.match_template(img_level, kern, pad_input=True, mode='constant', constant_values=0)
    local_maxi = peak_local_max(feat_img, indices=False, footprint=disk(4), exclude_border=False, threshold_rel=.4)
    
    # apply watershed
    markers = ndimage.label(local_maxi)[0]
    labels = watershed(-feat_img, markers, mask=img_gthres)
    
    # exclude regions that are too small, too large
    roi_props = skimage.measure.regionprops(labels, intensity_image=None, cache=True)
    bKeep = np.array([(rp['area']>5) & (rp['area']<70) & (rp['eccentricity']<.9) for rp in roi_props]).astype(bool)
    labels = np.zeros_like(labels)
    from itertools import compress
    for nLabel, rp in enumerate(compress(roi_props,bKeep)):
        labels[rp['coords'][:,0], rp['coords'][:,1]] = nLabel
    roi_props = skimage.measure.regionprops(labels, intensity_image=None, cache=True)

    # filter rois according to properities
    label_img = skimage.color.label2rgb(labels, image=None, colors=None, alpha=1, bg_label=0, bg_color=[0,0,0], image_alpha=.5)

    intermediates = {}
    intermediates['img_gamma'] = img_gamma
    intermediates['img_level'] = img_level
    intermediates['img_gthres'] = img_gthres
    intermediates['img_gthres'] = img_gthres
    rois = labels
    roi_image = label_img
    
    return roi_props, rois, roi_image, intermediates
"""
segment_nuclei = py"segment_nuclei"

##
segment_nuclei(img)

##
###############
## SEGMENTATION
###############
##
threshold = otsu_threshold(float.(h2b_zbrain))
##
imshow(h2b_zbrain.>threshold)
##
seg = fast_scanning(h2b_zbrain[:,:,100],threshold)
seg
##
