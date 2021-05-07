using ImageTracking, Images, ImageView, CoordinateTransformations, StaticArrays,
    LinearAlgebra, Glob, VideoIO, StatsBase, ProgressMeter, ImageSegmentation
    ImageFiltering, Random

video_dir = "/mnt/deissero/users/tyler/osprey/2021-04-21_lsd/pre-admin"
vid = openvideo(glob("*.mp4", video_dir)[1]);

##
function read_grayscale(vid)
    img = read(vid)
    channelview(img)[1,:,:]
end

img = read_grayscale(vid);
imshow(first_img)
H, W = size(img)
mean_frame = zeros(size(img)...);
##
seekstart(vid)
num_frames = 0
prog = ProgressUnknown("Frames read:")
while !eof(vid)
    img = read_grayscale(vid)
    num_frames += 1
    mean_frame .+= img
    ProgressMeter.next!(prog)
end;
ProgressMeter.finish!(prog)
##
mean_frame ./= num_frames
##
Gray.(mean_frame)
##
Gray.(img-mean_frame)
##
img_delta = img-mean_frame
bw = zeros(Bool,size(img))
bw .= false
bw[findall(img_delta .>= otsu_threshold(img_delta))] .= true
# bw[findall(img .>= 0)] .= 0
Gray.(bw)
##
erodeN(n) = n == 1 ? erode : erode ∘ erodeN(n-1)
dilateN(n) = n == 1 ? dilate : dilate ∘ dilateN(n-1)
bw = (~).(opening((~).(bw)))
Gray.(bw)
##
dist = 1 .- distance_transform(feature_transform(bw));
imshow(dist)
##
function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end
markers = label_components(dist .<= -0.5); sum(markers)
markers = label_components(bw .== 0);
##
segments = watershed(dist, markers) # can we speed this up...?5
imshow(map(i->get_random_color(i), labels_map(segments)) .* (1 .-bw)) # sanity check
##
# segments = unseeded_region_growing(Gray.(Float32.(bw)), 0.05);
