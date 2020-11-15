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
piezoDir = "/mnt/deissero/users/tyler/b115/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/ZSeries-piezo-2um-150um-045"
piezo2Dir = "/mnt/deissero/users/tyler/b115/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/ZSeries-piezo-2um-150um-second-046"
etlDir = "/mnt/deissero/users/tyler/b115/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/ZSeries-etl-2um-150um-044"

piezo = readZseriesTiffDir(piezoDir)
piezo2 = readZseriesTiffDir(piezo2Dir)
etl = readZseriesTiffDir(etlDir)
## WIP
# function julia_mutual_information(A, B)
#     "Mutual information for joint histogram"
#     # https://matthew-brett.github.io/teaching/mutual_information.html
#     h = fit(Histogram, (float(A), float(B)),nbins=20)
#     # Convert bins counts to probability values

#     pxy = h.weights / sum(h.weights)
#     px = sum(pxy, dims=2) # marginal for x over y
#     py = sum(pxy, dims=1) # marginal for y over x
#     px_py = reshape(px,:,1) * reshape(py,1,:)# Broadcast to multiply marginals
#     # Now we can do the calculation using the pxy, px_py 2D arrays
#     nzs = pxy .> 0 # Only non-zero pxy values contribute to the sum
#     return sum(pxy[nzs] .* log.(pxy[nzs] ./ px_py[nzs]))
# end

# julia_mutual_information(piezo[:,:,1], piezo2[:,:,1])
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

nZ = size(piezo,3)
mi = zeros(nZ, nZ)
miControl = zeros(nZ, nZ)

for pZ in 1:nZ
    for eZ in 1:nZ
        p = float(piezo[:,:,pZ])
        e = float(etl[:,:,eZ])
        c = float(piezo2[:,:,eZ])
        mi[pZ,eZ] = mutual_information(p, e)
        miControl[pZ,eZ] = mutual_information(p, c)
    end
end

##
m1 = copy(cartesianIndToArray(dropdims(argmax(miControl, dims=1), dims=1)))[:,1]
# m2 = copy(cartesianIndToArray(dropdims(argmax(mi, dims=2), dims=2)))[:,2]
# zs = collect(-50:50)
zs = collect(-50:2:100)
fig = plot(zs, 2 .* m1 .- 51, legend=:bottomright, label="")
plot!(zs, zs,label="")
title!("Piezo (2nd acquisition) planes with maximal MI to Piezo planes")
xlabel!("Piezo plane (μm)")
ylabel!("Piezo2 Plane (μm)")
savefig(fig, "figs/piezo_piezo_150um_max_mi.png")
fig
##
fig = heatmap(zs, zs, miControl)
xlabel!("Piezo plane (μm)")
ylabel!("Piezo2 Plane (μm)")
title!("Mutual information")
savefig(fig, "figs/piezo_piezo_150um_mi.png")
fig
##
m1 = copy(cartesianIndToArray(dropdims(argmax(mi, dims=1), dims=1)))[:,1]
# m2 = copy(cartesianIndToArray(dropdims(argmax(mi, dims=2), dims=2)))[:,2]
# zs = collect(-50:50)
zs = collect(-50:2:100)
fig = plot(zs, 2 .* m1 .- 51, legend=:bottomright, label="")
plot!(zs, zs,label="")
title!("Piezo planes with maximal MI to ETL planes")
xlabel!("ETL plane (μm)")
ylabel!("Piezo Plane (μm)")
savefig(fig, "figs/etl_piezo_150um_max_mi.png")
fig
##
fig = heatmap(zs, zs, mi)
xlabel!("ETL plane (μm)")
ylabel!("Piezo Plane (μm)")
title!("Mutual information")
savefig(fig, "figs/etl_piezo_150um_mi.png")
fig