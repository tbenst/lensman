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

"""Compute and assign local variables.
```
@experiment (etlVals, Z, H, W) = "2021-08-03/fish1/TSeries-1"
```

will expand to:
```
etlVals, Z, H, W = compute("2021-08-03/fish1/TSeries-1",
    ["etlVals", "Z", "H", "W"])
```
"""
macro experiment(ex)
    @assert ex.head == :(=) "no `=` found in expression."
    vars = ex.args[1] # :((etlVals, Z, H, W))
    exp_name = ex.args[2] # "2021-08-03/fish1/TSeries-1"
    vars_str = map(string, vars.args) # ["etlVals", "Z", "H", "W"]
    esc(:($vars = compute($exp_name, $vars_str)))
    # ex
