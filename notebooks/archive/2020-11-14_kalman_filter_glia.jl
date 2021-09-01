ENV["DISPLAY"] = "localhost:12"
using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕ # independent sum of Gaussian r.v.
using Statistics, Plots, Lensman, ImageView, Flux, CUDA
import Distributions, IterTools
using Parameters: @with_kw
using Flux.Data: DataLoader
using Base.Iterators: partition

using FileIO, NRRD, ImageView, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    ImageView, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, ImageView, Images, ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates

"ImageJ 'kalman filter'"
function kalmanFilter(ys::Array{<:Real,1};gain=0.8,useBurnIn=false)
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

tifdir = "/scratch/2020-10-06_gfap_h2b6s/fov1_hab/Tseries_resting-028"

glia = readZseriesTiffDir(tifdir, contains="Ch2");
neurons = readZseriesTiffDir(tifdir, contains="Ch3");
##
imshow(glia)

glia_float = Float64.(glia[1:10,1:10,1:1000])
##
# warmup JIT
filt = mapslices(kalmanFilter, glia[1:10,1:10,1:1000], dims=[3]);
filt = mapslices(kalmanFilter, neurons[:,:,1:1000], dims=[3]);
##

imshow(imadjustintensity(filt[:,:,10:1000]))

##
imshow(imadjustintensity(neurons[:,:,10:1000]))


## Flux LSTM
model = Chain(
    LSTM(1, 1),
    Dense(1,1,relu)
) # |> gpu

@with_kw mutable struct Args
    lr::Float64 = 3e-4
    bz::Int = 1024
    epochs::Int = 10
    seqlen::Int = 50
    device::Function = gpu
end

function images2batches(images; bz=64, seqlen=50, pad_val=0.)
    batches = []
    H, W = size(images)[1:2]
    for (i,j) in IterTools.product(1:H,1:W)
        # num_batches x bz x (seqlen)
        # bs = collect(partition(Flux.batchseq(Flux.chunk(images[i,j,:], bz),
        #     pad_val), seqlen))
        # bs = collect(partition(Flux.chunk(images[i,j,:], bz), seqlen))
        bs = collect(partition(Float32.(images[i,j,:]), seqlen))
        push!(batches, bs)
    end
    vcat(batches...)
end

# function train(images; kws...)
#     args = Args(; kws...)

#     train_data = 0
# end    
train_data = images2batches(glia[1:20,1:20,1:1000]);
# size(train_data), size(train_data[1]), size(train_data[1][1]), size(train_data[1][1][1]), size(train_data[1][1][1][1])

out = model.(train_data)
size(out), size(out[1])

## EXAMPLE
model2 = Chain(
  Dense(10, 5, σ),
  Dense(5, 2),
  softmax) |> gpu

batch_size = 5
model2.([rand(10) for i in 1:batch_size] |> gpu ) # => 5 x 2-element vector

## seq..?
seq = [rand(1) for i = 1:32] |> gpu
m = LSTM(1,2) |> gpu
out = m.(seq)
size(out), size(out[1]), size(out[1][1])

Flux.reset!(m)
@assert all(m(seq[2]) .!= out[2])
Flux.reset!(m)
@assert all(m(seq[1]) .== out[1])
@assert all(m(seq[2]) .== out[2])

## seq with batches...?
x_train = rand(1,50,128) # channel x time x nSamples
train_loader = DataLoader(x_train, batchsize=32) # channel x time x batch

m = LSTM(1,2)
@show m(x_train[:,1,1])
opt = ADAM()
loss(x,y) = Flux.Losses.mse(m(x),y)
nepochs = 5
Flux.train!(loss,Flux.params(m), IterTools.ncycle(train_loader, nepochs), opt)

## EXAMPLE
seq = [rand(10) for i = 1:32]
m = LSTM(10,2)
out = m.(seq)
size(out), size(out[1]), size(out[1][1])

## EXAMPLE
args = Args()


isfile("/tmp/input.txt") ||
        download("https://cs.stanford.edu/people/karpathy/char-rnn/shakespeare_input.txt","/tmp/input.txt")

text = collect(String(read("/tmp/input.txt")))

# an array of all unique characters
alphabet = [unique(text)..., '_']

text = map(ch -> Flux.onehot(ch, alphabet), text)
stop = Flux.onehot('_', alphabet)

N = length(alphabet)

# Partitioning the data as sequence of batches, which are then collected as array of batches
# num_batches x bz x (features)
Xs = collect(partition(Flux.batchseq(Flux.chunk(text, 25), stop), args.seqlen))
Ys = collect(partition(Flux.batchseq(Flux.chunk(text[2:end], 25), stop), args.seqlen))
size(Ys[1])