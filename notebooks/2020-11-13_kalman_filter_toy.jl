using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕ # independent sum of Gaussian r.v.
using Statistics, Plots
import Distributions

"ImageJ 'kalman filter'"
function kalmanFilter(ys;gain=0.8,useBurnIn=true)
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

"""
x[k] = Φx[k−1] + w[k],    w[k] ∼ N(0, Q)
y[k] = Hx[k] + v[k],    v[k] ∼ N(0, R)
"""
ts = collect(.1:.1:1000)
Gaussian(zeros(Int(size(ts,1)/4)),0.1)
xtrue = 1*sin.(ts);
# xtrue = cumsum(repeat(rand(Gaussian(zeros(Int(size(ts,1)/4)),0.1)), inner=4));
# xtrue = zeros(size(ts));
ys = rand(Gaussian(xtrue,0.5));
##
# priors
x0 = ys[1]

# we difference to stabilize but this doubles the variance
P0 = var(diff(ys))/2 # estimate the observation noise
diff1_var = var(diff(ys))
diff2_var = var(ys[3:end] - ys[1:end-2])
estimate = var(ys[3:end] - ys[1:end-2])/2 - P0

# estimate = (diff2_var - diff1_var)/2
actual = var(diff(xtrue))
actual2 = var(xtrue[3:end] - xtrue[1:end-2])

Φ = 1.
Q = 0.04 # state-change distribution
# @show Q = var(diff(xtrue)) # state-change distribution
H = 1.
R = P0



N = length(ys)
p = Gaussian(x0, P0)
ps = [p] # vector of filtered Gaussians
for i in 2:N
    global p
    # predict
    # p = Φ*p
    p = Φ*p ⊕ Gaussian(zero(x0), Q) #same as Gaussian(Φ*p.μ, Φ*p.Σ*Φ' + Q)
    # correct
    p, yres, _ = Kalman.correct(Kalman.JosephForm(), p, (Gaussian(ys[i], R), H))
    push!(ps, p) # save filtered density
end

l = @layout [a; b; c; d; e]
p1 = plot(ts, xtrue, legend=false, title="true")
p2 = plot(ts, ys, legend=false, title="obs")
kalmanPred = (x -> x.μ).(ps)
p3 = plot(ts, kalmanPred, legend=false, title="filtered")
imagejPred = kalmanFilter(ys,gain=.758)
p4 = plot(ts, imagejPred, legend=false, title="ImageJ 'kalman'")
p5 = plot(ts, kalmanFilter(ys, useBurnIn=false), legend=false, title="ImageJ (no burn in)")
# plot(p1, p2, p3, p4, p5, layout = l)
l = @layout [a; b; c;]
plot(p1, p3, p4, layout = l)
@show norm(xtrue-kalmanPred,2)
@show norm(xtrue-imagejPred,2)
