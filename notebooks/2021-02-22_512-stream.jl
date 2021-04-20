##
# ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose
import Gadfly
using Unitful: μm, m, s
## debug
using Sockets
HOSTNAME = "localhost"
PORT = 2001
pv = connect(HOSTNAME, PORT)
println("connected")
data = []

@async while isopen(pv)
    println("open")
    push!(samples, read(pv))
    println("pushed")
end

##
# HOSTNAME = "171.65.101.131" # B113
HOSTNAME = "171.65.101.179" # B115
PORT = 1050
pv = connect(HOSTNAME, PORT)
println("connected")

# First, we read sizes of images to follow
bSize = read(pv, 8)
# reshape fails unless an Int64
C, Z, H, W = Int64.(reinterpret(Int16,bSize))
# TODO: right now only z=1 is functional...???
# at z=4, I get >50% more samples than expected...?
# perhaps "bad" pixels arent' being invalidated...?
# -> but no consistent visual effect if so...
# Z = 3
@show C, Z, H, W


println("writing to host")
write(pv, "Tseries")

samples = []
n = 500000
@async while isopen(pv)
    bArray = read(pv, sizeof(Int16)*n)
    arr = reinterpret(Int16, bArray)
    push!(samples, arr)
end

## 1-plane 1024
expected_Z = 1
expected_W = 1024
expected_H = 1024
expected_C = 2
expected_rate = 15
expected_time = 30 # sec
# num_reps = expected_rate * expected_time
num_reps = 451
samples_per_rep = expected_W * expected_H * expected_Z * expected_C
expected = samples_per_rep * num_reps
actual = length(samples) * n
@show expected_Z, expected_W, expected_H, expected_C, expected_rate, expected_time, expected, actual 
print("fraction of expected samples: $(actual/expected)")
print("expected reps: $num_reps, actual reps: $(actual/samples_per_rep)")
# 1.004

## 10-plane 1024
expected_Z = 10
expected_W = 1024
expected_H = 1024
expected_C = 2
expected_rate = 1/0.666666
expected_time = 30 # sec
samples_per_rep = expected_W * expected_H * expected_Z * expected_C
# num_reps = expected_rate * expected_time
num_reps = 46
expected = samples_per_rep * num_reps
actual = length(samples) * n
@show expected_Z, expected_W, expected_H, expected_C, expected_rate, expected_time, expected, actual 
print("fraction of expected samples: $(actual/expected)")
print("expected reps: $num_reps, actual reps: $(actual/samples_per_rep)")
# 1.004


## 1-plane 512x512
expected_Z = 1
expected_W = 512
expected_H = 512
expected_C = 2
expected_rate = 30
expected_time = 30 # sec
samples_per_rep = expected_W * expected_H * expected_Z * expected_C
# num_reps = expected_rate * expected_time
num_reps = 89
expected = samples_per_rep * num_reps
actual = length(samples) * n
@show expected_Z, expected_W, expected_H, expected_C, expected_rate, expected_time, expected, actual 
print("fraction of expected samples: $(actual/expected)")
print("expected reps: $num_reps, actual reps: $(actual/samples_per_rep)")

## 10-plane 512x512
expected_Z = 10
expected_W = 512
expected_H = 512
expected_C = 2
expected_rate = 0.3333
expected_time = 30 # sec
samples_per_rep = expected_W * expected_H * expected_Z * expected_C
# num_reps = expected_rate * expected_time
num_reps = 89
expected = samples_per_rep * num_reps
actual = length(samples) * n
@show expected_Z, expected_W, expected_H, expected_C, expected_rate, expected_time, expected, actual 
print("fraction of expected samples: $(actual/expected)")
print("expected reps: $num_reps, actual reps: $(actual/samples_per_rep)")

## 1-plane 256x256
expected_Z = 10
expected_W = 256
expected_H = 256
expected_C = 2
expected_rate = 1/0.17
expected_time = 30 # sec
samples_per_rep = expected_W * expected_H * expected_Z * expected_C
# num_reps = expected_rate * expected_time
num_reps = 174
expected = samples_per_rep * num_reps
actual = length(samples) * n
@show expected_Z, expected_W, expected_H, expected_C, expected_rate, expected_time, expected, actual 
print("fraction of expected samples: $(actual/expected)")
print("expected reps: $num_reps, actual reps: $(actual/samples_per_rep)")

## 10-plane
# expected_Z = 10
# expected_W = 512
# expected_H = 512
# expected_C = 2
# expected_rate = 3
# expected_time = 30 # sec
# expected = expected_W * expected_H * expected_Z * expected_C * expected_rate * expected_time
# actual = length(samples) * n
# print("fraction of expected samples: $(actual/expected)")
##




##  10-plane reshaping test
array1d = vcat(samples...)[1:expected]
expected / (C * H * W * Z)
ims = reshape(array1d, C, W, H, Z, num_reps)
# ims = permutedims(ims, [1,2,1,4,5])
# ims = mean(ims, dims=5)
imshow(ims[:,:,1,:,:])


##

currentVolume = Observable(zeros(Int16,C,H,W,Z))

numFrames = Observable(0)

images = []

# callback for new volumes
on(currentVolume) do vol
    numFrames[] += 1
    println("got volume $(numFrames[])")
    push!(images, vol)
end

# on(currentVolume) do vol
#     println("got volume with mean $(mean(vol))!")
# end

# start acquisition by sending TCP
tmpVol = zeros(Int16,C,H,W,Z)
println("writing to host")
write(pv, "Tseries")
# first image is from acquisition prior to now...
_ = read(pv, sizeof(Int16)*C*H*W)
# Next, we start reading images as they come in...
# TODO we're dropping frames, why?? PV python server
# says "acquired frame 1013"
# and we received 1012, but supposed to get 1024...
@async while isopen(pv)
    for z in 1:Z
        bArray = read(pv, sizeof(Int16)*C*H*W)
        println("read z $z")
        arr = permutedims(reshape(reinterpret(Int16, bArray),C,W,H),[1,3,2])
        # we have to reassign currentVolume to trigger observable listener...
        tmpVol[:,:,:,z] .= arr
    end
    currentVolume[] = copy(tmpVol)
end

nothing
##
close(pv) # warning: only run once done!
function zeroAdjust(image::Array{T,4}) where T<:Real
    # PMT offset
    # @assert sum(image.==0) == 0
    im = image .- 8192
    im[im.<0] .= 0
    # every other line, reverse column for roundtrip correction
    # for B115, need to start with odd line else fliplr
    im[:,1:2:end,:,:] .= im[:,1:2:end,end:-1:1,:]
    im
end

# compile
zeroAdjust(images[1])

avgImage = mean(zeroAdjust.(images))[2,:,:,1]
imshow(adjust_gamma(avgImage,0.5))
