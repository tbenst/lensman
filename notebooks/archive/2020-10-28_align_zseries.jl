## Take two zseries at different wavelengths, and find planewise correspondance 
using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5
    import Base.Threads.@threads
##
# 100 um (+/- 50)
# piezoDir = "/mnt/deissero/users/tyler/b115/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/ZSeries-piezo-1um-037"
# etlDir = "/mnt/deissero/users/tyler/b115/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/ZSeries-etl-1um-038"

# 150um (-50 to 100)
w1Dir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/ZSeries-structural-049"
w2Dir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/ZSeries-structural-920-050"

zseries1 = readZseriesTiffDir(w1Dir)
zseries2 = readZseriesTiffDir(w2Dir)
##
py"""
import numpy as np
def mutual_information(A, B):
    "Mutual information for joint histogram"
    # https://matthew-brett.github.io/teaching/mutual_information.html
    hgram, x_edges, y_edges = np.histogram2d(A.ravel(), B.ravel(), bins=20) 
    # Convert bins counts to probability values

    pxy = hgram / float(np.sum(hgram))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
    # Now we can do the calculation using the pxy, px_py 2D arrays
    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))
"""
mutual_information = py"mutual_information"

nZ1 = size(zseries1,3)
nZ2 = size(zseries2,3)
mi = zeros(nZ1, nZ2)

for z1 in 1:nZ1
    for z2 in 1:nZ2
        p = float(zseries1[:,:,z1])
        q = float(zseries2[:,:,z2])
        mi[z1,z2] = mutual_information(p, q)
    end
end

##
m1 = copy(cartesianIndToArray(dropdims(argmax(mi, dims=1), dims=1)))[:,1]
zs = collect(1:nZ1)
fig = plot(zs, m1, legend=:bottomright, label="")
plot!(zs, zs,label="")
title!("840nm planes with maximal MI to 920nm plane")
xlabel!("plane 920nm")
ylabel!("plane 840nm")
savefig(fig, "figs/840nm_920nm_max_mi.png")
fig
##
fig = heatmap(zs, zs, mi)
xlabel!("plane 920nm")
ylabel!("plane 840nm")
title!("Mutual information")
savefig(fig, "figs/840nm_920nm_mi.png")
fig
