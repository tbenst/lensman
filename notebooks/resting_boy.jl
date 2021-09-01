##
ENV["DISPLAY"] = "localhost:10.0"
using ImageView, Images

# These could still be interesting because there was no movement between
# btf_path = "/data/allan/olympus/2021-02-17_elavl3-GCaMP6f_7dpf/test_fish2/tseries_1plane_20min_post_embedding.ome.btf"  # 14
btf_path = "/data/allan/olympus/2021-02-17_elavl3-GCaMP6f_7dpf/test_fish2/tseries_1plane_60min_post_embedding.ome.btf"  # 15

# Hopefully these are without movement
# btf_path = "/data/allan/olympus/2021-02-17_elavl3-GCaMP6f_7dpf/test_fish2/tseries_1plane_100min_post_embedding.ome.btf"  # 10
# btf_path = "/data/allan/olympus/2021-02-17_elavl3-GCaMP6f_7dpf/test_fish2/tseries_1plane_140min_post_embedding.ome.btf"  # 11
# btf_path = "/data/allan/olympus/2021-02-17_elavl3-GCaMP6f_7dpf/test_fish2/tseries_1plane_180min_post_embedding.ome.btf"  # 12
# btf_path = "/data/allan/olympus/2021-02-17_elavl3-GCaMP6f_7dpf/test_fish2/tseries_1plane_220min_post_embedding.ome.btf"  # 13

using PyCall
tifffile = pyimport("tifffile")
x = tifffile.imread(btf_path)  # (T, H, W)

cv2 = pyimport("cv2")
mask_path = "/home/allan/workspace/lensman/tseries_1plane_20min_post_embedding_lh.png"
# mask_path = "/home/allan/workspace/lensman/tseries_1plane_100min_post_embedding_lh.png"
mask = cv2.imread(mask_path)[1:end, 1:end, 2]
mask = reinterpret(Bool, mask)
mask = reshape(mask, (size(mask)..., 1))
mask = convert(Array{Bool,3}, mask)

using Lensman

(T, H, W) = size(x)
x = permutedims(x, [2, 3, 1])
x = reshape(x, (H, W, 1, T))

corr_image = imageCorrWithMask(x, mask)
println("$(minimum(corr_image)) $(maximum(corr_image))")

##
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
matplotlib.use("tkagg")
fig = plt.figure()
cmax = 0.3
plt.imshow(corr_image, cmap="RdBu_r", clim=(-cmax, cmax))
plt.colorbar()
# plt.show()
display(fig)

## Looking reasonable (make sure code correct), now do for other times

# This uses ImageMagick and seems to do soe messed up RGB thing
# ```
# x = load(btf_path)
# size(x)
# julia> typeof(x)
# Array{RGB{Normed{UInt16,16}},3}
# ```

# tif0 =  reinterpret(UInt16, ImageMagick.load(btf_path))
# imshow(imadjustintensity(tif0[1:end, 1:end, 500]))
