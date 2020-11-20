##
using Sockets, Observables, Statistics, Images, ImageView
``
ENV["DISPLAY"] = "localhost:12.0"
HOSTNAME = "171.65.101.131" # B113
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
Z = 3
@show C, Z, H, W

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
write(pv, "Tseries")
# first image is from acquisition prior to now...
_ = read(pv, sizeof(Int16)*C*H*W)
# Next, we start reading images as they come in...
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
function zeroAdjust(image::Array{<:Real,4})
    # PMT offset
    @assert sum(image.==0) == 0
    # im = image .- 8192
    # im[im.<0] .= 0
    # every other line, reverse column for roundtrip correction
    im[:,2:2:end,:,:] .= im[:,2:2:end,end:-1:1,:]
    im
end
imshow(images[3][2,:,:,1])
# zIm = zeroAdjust.(images);
# imshow(adjust_gamma(zIm[1][2,:,:,:],0.5))
# imshow(adjust_gamma(mean(zIm)[2,:,:,:],0.5))
##

# vol = reshape(currentVolume[], H,W,Z,C)
imshow(currentVolume[][1,:,:,:])
println(numFrames[])

##
imshow(rand(16,16,4,2))