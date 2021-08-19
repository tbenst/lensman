"Save .h5 file with left & right hemishpere masks."
using AlgebraOfGraphics, CairoMakie
using Lensman, Images, Glob, NPZ, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      H5Sparse, SparseArrays, CategoricalArrays
import Lensman: @pun, @assign
using ImageView
# import PyPlot
# plt = PyPlot
# matplotlib = plt.matplotlib
# np = pyimport("numpy")
# import Plots
# import Plots: heatmap
# using Gadfly
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
aog = AlgebraOfGraphics
Data = aog.data
set_aog_theme!() #src
##
# resources = Resources();
##
@pun (h2b_zbrain, zbrain_masks, zbrain_mask_names, zbrain_dir) = resources
# imshow(h2b_zbrain)
##
hprofile = mean(collect(h2b_zbrain), dims=(2,3))[:,1,1];
data((x=hprofile,)) * mapping(:x) * visual(Lines) |> draw

##
avgcols = mean(collect(h2b_zbrain), dims=2)[:,1,:]
avgdepth = mean(collect(h2b_zbrain), dims=3)[:,:,1]
##
@assert size(avgcols) == size(h2b_zbrain)[[1,3]] # H, Z
optH_byZ =  mapslices(argmax, avgcols, dims=1)[1,:];
data((x=optH_byZ,)) * mapping(:x) * visual(Lines) |> draw
##
"Convert matrix (wide data) to tidy (tall data)."
function tidy_matrix(matrix)
    matrix_df = DataFrame(matrix, string.(collect(1:size(matrix,2))))
    matrix_df[!,:y] = 1:size(matrix,1)
    matrix_df = stack(matrix_df, 1:size(matrix,2))
    rename!(matrix_df, :variable => :x)
    matrix_df[!, :x] = parse.(Int,matrix_df[!, :x])
    matrix_df
end

function tidy_3d(tensor)
    df = DataFrame()
    ys = Int[]
    xs = Int[]
    zs = Int[]
    values = Float64[]
    for i in 1:size(tensor,1)
        for j in 1:size(tensor,2)
            for k in 1:size(tensor,3)
                push!(ys, i)
                push!(xs, j)
                push!(zs, k)
                push!(values, tensor[i,j,k])
            end
        end
    end
    DataFrame(y=ys, x=xs, z=zs, value=values)
end
tidy_3d(h2b_zbrain[1:10,1:10,1:10])
h2b_df = tidy_3d(h2b_zbrain)

##
avgdepth_df = tidy_matrix(avgdepth)
avgcols_df = tidy_matrix(avgcols)
# df = tidy_matrix(avgcols)
# df2 = (x=[0, 138], y=[0.5, 0.5])
# df = tidy_matrix(1 ./ (1 .+ avgcols[311-50:311+50,:])) 
# df[:,:value] = log.(df[:,:value].+1)

# p = data(df)  * mapping(:x, :y, :value) * (
#      visual(Heatmap)
# ) + data(df2) * visual(Lines)|> draw

using Polynomials, LsqFit
using FixedPolynomials
import DynamicPolynomials: @polyvar, polynomial, monomials
##
# model(x,p) = Polynomial(p).(x)
# p0 = [311, 0.5, 0.5, 0.5]
# p0 = [300, 0.03, -7.8e-5, 4.4e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
function make_mvar_poly(deg, coeff)
    @polyvar x y
    vars = monomials([x,y], 1:deg)
    polynomial(coeff[1:end-1], vars) + coeff[end]
end

fast_deg3(a...) = (x,y)->a[1]+a[2]*x+a[3]*y+a[4]*x*y+a[5]*x^2+a[6]*y^2+a[7]*x^2*y+a[8]*x*y^2+a[9]*x^3+a[10]*y^3
# f = make_mvar_poly(3, collect(1:10))
f = fast_deg3(collect(1:10)...)
f(2,3)
##
# model(x,p) = make_mvar_poly(3, p).(x[:,1], x[:,2])
model(x,p) = fast_deg3(p...).(x[:,1], x[:,2])
p0 = [310, (randn(9) .* map(x->10.0^(-x), 1:9))...]
# fit = curve_fit(model, df[:,:x], df[:,:y], df[:,:value], p0)
##
# df = h2b_df[1:10000:end,:]
# df = h2b_df[1:100:end,:]
df = h2b_df
fit = curve_fit(model, Array(df[:,[:x,:z]]), df[:,:y], df[:,:value], p0)
# fit = curve_fit(model, Array(df[:,[:x,:z]]), df[:,:y], df[:,:value], p0)

# xs = collect
# best = Polynomial(fit.param)
fit.param
# deg 3 x,y:
# 313.8306435446811
# 0.03162144712562135
# -0.07162310094386597
# 8.769285377672171e-5
# -7.588491273339297e-5
# -0.001414943301453414
# 2.1580350392756406e-7
# -2.6683579146412668e-6
# 3.765100092436289e-8
# 1.583931520361477e-5
##
best = fast_deg3(fit.param...)
xs = collect(1:size(h2b_zbrain,2))
zs = collect(1:size(h2b_zbrain,3))

##
# using `best` function, plot the line for each 
# save as a boolean mask
centerline_mask = zeros(Bool, size(h2b_zbrain)...)
for x in 1:size(h2b_zbrain,2)
    for z in 1:size(h2b_zbrain,3)
        y = Int(round(best(xs[x],zs[z])))
        centerline_mask[y,x,z] = true
    end
end
im = RGB.(h2b_zbrain)
selectdim(channelview(im),1,1) .= centerline_mask
imshow(im)

##
"""
Calculate running intersection of all masks.
"""
function intersect_masks(zbrain_masks)
    H = Int(zbrain_masks["height"])
    W = Int(zbrain_masks["width"])
    Z = Int(zbrain_masks["Zs"])
    mask = zeros(Bool, H, W, Z)
    @showprogress for j in 1:size(zbrain_masks["MaskDatabase"],2)
        new = reshape(Array(zbrain_masks["MaskDatabase"][:,j]), (H, W, Z))
        # new = permutedims(new, (2, 1, 3))
        mask = mask .| new
    end
    mask
end
brain_mask = intersect_masks(zbrain_masks)
##
left_mask = zeros(Bool, size(h2b_zbrain)...)
right_mask = zeros(Bool, size(h2b_zbrain)...)
ys = collect(1:size(h2b_zbrain,1))
for x in 1:size(h2b_zbrain,2)
    for z in 1:size(h2b_zbrain,3)
        y = Int(round(best(xs[x],zs[z])))
        # by manual inspection, left is a little shortchanged
        left_mask[:,x,z] .= ys .<= y
        right_mask[:,x,z] .= ys .> y
    end
end
revert_to_zbrain_orientation(m) = permutedims(reverse(reverse(m,dims=2),dims=3), (2, 1, 3))

left_mask = revert_to_zbrain_orientation(left_mask) .& brain_mask
right_mask = revert_to_zbrain_orientation(right_mask) .& brain_mask;
##
im2 = RGB.(imadjustintensity(revert_to_zbrain_orientation(h2b_zbrain)))
# im2 = RGB.(left_mask)
channelview(im2)[1,:,:,:] .= left_mask
channelview(im2)[3,:,:,:] .= right_mask
imshow(im2)

##
using H5Sparse
h5path = "$zbrain_dir/hemisphere_masks.h5"
h5 = h5open(h5path, "w")
##
# h5["size"] = collect(size(left_mask))
mask2vec(x) = sparse(reshape(x, prod(size(left_mask)),1))
H5SparseMatrixCSC(h5, "left", mask2vec(left_mask))
H5SparseMatrixCSC(h5, "right", mask2vec(right_mask))
close(h5)