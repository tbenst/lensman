# code that is not used but I can't bear parting with

function adj_gamma_rgb(img; gamma=0.5)
    colorview(RGB, adjust_histogram(channelview(img), GammaCorrection(gamma)))
end

"combine three seperate 1-channel RGB into a 3-channel RGB."
function singlechanrgb2rgb(red, green, blue)
    new = copy(red)
    greenview = selectdim(channelview(new),1,2)
    greenview[:] = selectdim(channelview(green),1,2)
    blueview = selectdim(channelview(new),1,3)
    blueview[:] = selectdim(channelview(blue),1,3)
    new
end

"layer ontop of base."
function layer_imgs(base, layer; ratio=0.8)
    new = base .* (1-ratio)
    new += layer .* ratio
    new
end
