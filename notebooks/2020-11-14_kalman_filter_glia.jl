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
