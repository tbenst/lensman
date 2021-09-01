using Random

function get_random_color(seed::Float64)
    seed = reinterpret(Int64, seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
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