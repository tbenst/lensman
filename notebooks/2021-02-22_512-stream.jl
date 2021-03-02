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

push!(samples, read(pv))
@async while isopen(pv)
    push!(samples, read(pv))
end
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
        # we have to reassign currentVolume to trigger observable listener...
        arr = permutedims(reshape(reinterpret(Int16, bArray),C,W,H),[1,3,2])
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
