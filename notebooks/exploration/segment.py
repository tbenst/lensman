#%%
import os
import numpy as np
import scipy
from scipy import ndimage
import skimage
import skimage.io
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


datadir = "/data/dlab/zfish_bruker/2020-08-04_elavl3-ChRmine-Kv2.1_elavl3-H2B-GCaMP6s/fish1/fov1-hab-brainwide/green-003/"
datadir = r"E:\tyler\2020-10-12_elavl3-chrmine-kv2.1_h2b-6s_11dpf\fov2_hab_raphe\SingleImage-820_128avg-031"
imgpath = glob(os.path.join(datadir,"green*.tif"))[0]
img = skimage.io.imread(imgpath)

seg = segment_nuclei(img)
print("yes")
# %%
seg
# %%
from glob import glob
import matplotlib.pyplot as plt
plt.imshow(seg[2])
# %%
