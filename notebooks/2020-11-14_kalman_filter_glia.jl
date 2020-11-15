using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕ # independent sum of Gaussian r.v.
using Statistics, Plots
import Distributions

"ImageJ 'kalman filter'"
function kalmanFilter(ys;gain=0.8,useBurnIn=false)
    predictions = zeros(size(ys))
    predictions[1] = ys[1]
    for t in 1:size(ys,1)-1
        K = useBurnIn ? 1 / (t+1) : 0.
        pred = predictions[t]
        obs = ys[t+1]
        predictions[t+1] = pred*gain + obs*(1-gain) + K*(obs-pred)
    end
    predictions
end

tifdir = "/mnt/deissero/users/tyler/b113/2020-11-04_elavl3-chrmine-Kv2.1_h2b6s_7dpf/SingleImage-840nm-010"

fishDir, _ = splitdir(tifdir)
tiff_files = joinpath.(tifdir, filter(x->(x[end-6:end]=="ome.tif") & occursin("Ch3", x),
    readdir(tifdir)))
plane = ImageMagick.load(tiff_files[1])
