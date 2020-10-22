# TAKEAWAY: move 820nm down by 2μm to equal 920nm for 25X on Bruker B115

##
using Plots, ImageMagick, Lensman, Images
##
dataDir = "/mnt/deissero/users/tyler/b115/20201020_RegisterWavelengths"
dir820, _, _, dir920 = readdir(dataDir)

get_green_images(x) = (x[end-3:end]==".tif") & occursin("Ch3", x)

tifdir820 = joinpath(dataDir,dir820)
tifpaths820 = joinpath.(tifdir820,filter(get_green_images,
    readdir(tifdir820)))
stack820 = read_zseries(tifpaths820)

tifdir920 = joinpath(dataDir,dir920)
tifpaths920 = joinpath.(tifdir920,filter(get_green_images,
    readdir(tifdir920)))
stack920 = read_zseries(tifpaths920);


##
N = size(stack920,3)

normed920 = imadjustintensity(sum(stack920, dims=(1,2))[1,1,:])
plot(collect(1:N), normed920, label="920")

normed820 = imadjustintensity(sum(stack820, dims=(1,2))[1,1,:])
plot!(collect(1:N), normed820, label="820")

maxshift = 10
for shift in 1:maxshift
    println("lshift $shift: $(transpose(normed820[1+shift:end-maxshift+shift])*normed920[1+maxshift-shift:end-shift])")
end
for shift in 1:maxshift
    println("rshift $shift: $(transpose(normed920[1+shift:end-maxshift+shift])*normed820[1+maxshift-shift:end-shift])")
end

##
shift = 2
normed920 = imadjustintensity(sum(stack920, dims=(1,2))[1,1,:])
plot(collect(1:N-shift), normed920[1+shift:end], label="920")

normed820 = imadjustintensity(sum(stack820, dims=(1,2))[1,1,:])
plot!(collect(1:N-shift), normed820[1:end-shift], label="820")
