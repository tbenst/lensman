using Lensman
using Statistics, ImageFiltering, JuMP
L = Lensman
raw_path = L.get_series_dir("2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123",
 L.DEFAULT_SETTINGS[:tseries_root_dirs])
deinterlaced_path = L.get_series_dir("2021-06-08_rsChRmine_h2b6s/fish2/deinterlacing/TSeries-lrhab-titration-123/Deinterlaced",
 L.DEFAULT_SETTINGS[:tseries_root_dirs])

# Prairie view block ripper deinterlacing is almost worthless
raw = avgTseries(raw_path);
deinterlaced = avgTseries(deinterlaced_path);
Gray.(hcat(
    imadjustintensity(raw[:,:,6]),
    imadjustintensity(deinterlaced[:,:,6])
))

##

function quick_show(img)
    # Gray.(adjust_histogram(imadjustintensity(img), Equalization()))
    Gray.(adjust_histogram(imadjustintensity(img), GammaCorrection(0.2)))
end
"""Average vertical variance of 3 rows for image stack.

inner means skip top & bottom (e.g. outside of brain)."""
function interlace_variance(img; inner=0)
    n = ndims(img)
    window = ones(Int,n)
    window[2] = 3
    @show inner:(lastindex(img,2)-inner)
    inner = selectdim(img,2,(inner+1):(lastindex(img,2)-inner))
    mean(mapwindow(var, img, Tuple(window)))
end

"function is worthless...?"
function combing_metric(vec)
    small = abs.(diff(x))
    big = vcat(abs.(diff(x[1:2:end])), abs.(diff(x[2:2:end])))
    small = sum(small)/length(small)
    big = sum(big)/length(big)
    combing = big - small
end

"""Average vertical variance of 3 rows for image stack.

inner means skip top & bottom (e.g. outside of brain)."""
function interlace_difference(img; inner=0)
    n = ndims(img)
    # window = ones(Int,n)
    # window[2] = 3
    # window[2] = 5
    # inner = selectdim(img,1,(inner+1):(lastindex(img,1)-inner))
    # f = x->abs(x[2]-x[1])
    # expect neighbor to be large and two off to be small
    # if combing this is large, no combing closer to 0
    # mean(mapwindow(combing_metric, img, Tuple(window)))
    # every2 = mean(abs.(diff(selectdim(raw,2,1:2:lastindex(raw,2)),dims=2)))
    every1 = mean(abs.(diff(raw,dims=1)))
    # every2 - every1
    every1
end

"Shift vector and zero pad"
function shift_vector(vector, shift::Int)
    maxidx = length(vector)
    the_shift = maximum([1, 1+shift]):minimum([maxidx, maxidx+shift])
    new = zero(vector)
    if shift >= 0
        new[the_shift] = vector[1:length(the_shift)]
    else
        new[the_shift] = vector[(maxidx-length(the_shift)+1):end]
    end
    new
end

x = collect(1:10)
@assert shift_vector(x,3) == [0,0,0,1,2,3,4,5,6,7]
@assert shift_vector(x,-3) == [4,5,6,7,8,9,10,0,0,0]
@assert shift_vector(x,0) == x

"Shift even rows by `shift` pixels for nd image."
function shift_every_other_row(img,shift::Integer; circular=false, inner=0)
    inner_img = selectdim(img,2,(inner+1):(lastindex(img,2)-inner))
    new = zero(inner_img)
    # select odd rows
    new_odds = selectdim(new,1,1:2:lastindex(inner_img,1))
    # assign unchanged
    new_odds .= selectdim(inner_img,1,1:2:lastindex(inner_img,1))
    
    new_evens = selectdim(new,1,2:2:lastindex(inner_img,1))
    # shift even rows
    evens = selectdim(inner_img,1,2:2:lastindex(inner_img,1))
    if circular
        f = vec -> circshift(vec, shift)
    else
        f = vec -> shift_vector(vec,shift)
    end
    shifted_evens = mapslices(f, evens,dims=2)
    new_evens .= shifted_evens
    new
end
# quick_show(shift_every_other_row(raw[:,:,6],-2))

# @show interlace_variance(shift_every_other_row(raw,-2))
# @show interlace_variance(raw)
shifts = collect(-5:5)
# scores = map(x->interlace_variance(
scores = map(x->interlace_difference(
    shift_every_other_row(raw,x, circular=false,inner=0)),
    shifts)
println("best shift: $(shifts[argmin(scores)])")
scores
# model = Model(GLPK.Optimizer)
##
# manual looks good, but automated methods failing
quick_show(shift_every_other_row(raw,-2, circular=false,inner=0))[:,:,6]

interlace_difference(shift_every_other_row(raw,0, circular=false,inner=0))