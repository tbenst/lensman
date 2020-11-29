##
ENV["DISPLAY"] = "localhost:10.0"
using Sockets, Observables, Statistics, Images, ImageView, Lensman

HOSTNAME = "171.65.101.131" # B115
# HOSTNAME = "171.65.101.179" # B115
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
    im = image .- 8192
    im[im.<0] .= 0
    # every other line, reverse column for roundtrip correction
    # for B115, need to start with odd line else fliplr
    im[:,1:2:end,:,:] .= im[:,1:2:end,end:-1:1,:]
    im
end
avgImage = mean(zeroAdjust.(images))[2,:,:,1]
imshow(adjust_gamma(avgImage,0.5))
##
candNeuronCenterMask = findNeurons(avgImage,1.3, 15);
Gray.(candNeuronCenterMask)
##
candidateTargetLocs = findall(candNeuronCenterMask)
candidateTargets = copy(mask2ind(candNeuronCenterMask))

##
img = RGB.(imadjustintensity(avgImage))
stim_points = zeros(Bool,size(img))
neuron_locs = candidateTargets
stim_points[neuron_locs] .= true
stim_points = dilate(dilate(stim_points))
channelview(img)[[1,3],:,:,:] .= 0
channelview(img)[1,:,:,:] .= 0.5*float(stim_points)
img


##
using Lensman
using Unitful: μm, m, s, uconvert

## to burn at etl=0 if using calibration circa 2020-11-16, need +30 offset
offset = float(uconvert(m, 30μm)) / m * 2
targets1 = [280. 553. offset;
            661. 486. offset;
            642. 422. offset;
            715. 477. offset;] ./ 2

create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-24_quick_h2b6f")
