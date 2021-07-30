# MUST RUN ON LINUX FOR ANTS support (unless compiled for windows...)

# This notebook takes a 820nm Zseries, and:
# - registers to Zbrain
# - finds the imaging planes correspondence to the z-series
##
# ENV["DISPLAY"] = "localhost:12.0"
using FileIO, NRRD, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations, TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, Images,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, Lensman, ImageFiltering, Glob, Plots, HDF5,
    Dates, ProgressMeter, CSV, DataFrames
import Base.Threads.@threads
using Unitful: μm, m, s, mW
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib

# using ImageView
using Unitful: μm, m, s
sio = pyimport("scipy.io")

if Sys.iswindows()
    @error "this script only supports linux!"
end

ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
    zbrain_dir = "/oak/stanford/groups/deissero/users/tyler/zbrain"
    slmDir = "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files"
else
    # tseriesRootDir = "/data/dlab/b115"
    tseriesRootDir = "/scratch/b115"
    zbrain_dir = "/mnt/deissero/users/tyler/zbrain"
    ants_path = "/opt/ANTs/install/bin/antsRegistration"
    slmDir = "/mnt/deissero/users/tyler/b115/SLM_files"
    # tseriesRootDir = "/mnt/deissero/users/tyler/b115"
end

##
zseriesdir = joinpath(tseriesRootDir, "2021-07-14_rsChRmine_h2b6s_5dpf/fish1",
    "ZSeries-structural-840nm-048")

tseriesdir = joinpath(tseriesRootDir, "2021-07-14_rsChRmine_h2b6s_5dpf/fish1",
    "TSeries-lrhab-118trial-061")

fishDir, _ = splitdir(zseriesdir)

zbrain_units = (0.798μm, 0.798μm, 2μm)

# Read Zseries from file
# zseriesdir = joinpath(pv.directory,"Zseries_024034.17_-000")
tiff_files = joinpath.(zseriesdir, filter(x->(x[end-6:end]=="ome.tif") & occursin("Ch3", x),
    readdir(zseriesdir)))
tif0 = ImageMagick.load(tiff_files[1])
zseriesH, zseriesW = size(tif0)

zseries = zeros(Normed{UInt16,16}, zseriesH, zseriesW, size(tiff_files, 1))
@threads for z in 1:size(tiff_files,1)
    zseries[:,:,z] = ImageMagick.load(tiff_files[z])
end
size(zseries)

#####
    ## Read microns Per pixel (from file...)
#####

# open(glob("*.xml",zseriesdir)[1], "r") do io
#     zseries_xml = read(io, String)
#     global zseries_xml = xp_parse(zseries_xml)
# end

# "apply each function to args"
# funprod(functions...) = (args...)->map(x->x(args...), functions)

# micronsPerPixel_xml = zseries_xml[xpath"""//PVStateValue[@key="micronsPerPixel"]"""][1]
# # parse xml
# lookup_μm = axis -> etree -> parse(Float64,
#     micronsPerPixel_xml[xpath"""IndexedValue[@index="$axis"]"""][1].attr["value"])μm
# microscope_units = funprod(map(lookup_μm, ["YAxis", "XAxis", "ZAxis"])...)(micronsPerPixel_xml)

# convert to 25x units (16x is actually 14.4x, assumes Prairie View has 16x selected)
# @info "assuming 16x units but using 25x objective"
# microscope_units = ((microscope_units[1:2] .* (14.4/25))..., microscope_units[3])
# microscope_units = (0.5μm, 0.5μm, 1.5μm)

# TODO: make Olympus 25x in bruker have proper units....
@info "hardcoded z μ unit"
zseries_lateral_unit = microscope_lateral_unit(zseriesW)
zseries_microscope_units = (zseries_lateral_unit, zseries_lateral_unit, 2.0μm)

# zseries = centered(zseries) # need centered for qd registration
zseries = AxisArray(zseries, (:y, :x, :z), zseries_microscope_units)
adj_zseries = adjust_histogram(imadjustintensity(zseries), Equalization())
adj_zseries = AxisArray(adj_zseries, (:y, :x, :z), zseries_microscope_units)
@assert zseries_microscope_units[1] > 0.1μm
zseries_microscope_units

#####
# Read Zbrain H2B
#####
h2b_zbrain = AxisArray(permutedims(
        h5read("$zbrain_dir/AnatomyLabelDatabase.hdf5", "Elavl3-H2BRFP_6dpf_MeanImageOf10Fish"),
        (2,1,3)),
    (:y, :x, :z), zbrain_units)
ps = pixelspacing(h2b_zbrain); SD = ps./maximum(ps)
adj_h2b_zbrain = adjust_histogram(imadjustintensity(h2b_zbrain), Equalization())
adj_h2b_zbrain = AxisArray(adj_h2b_zbrain, (:y, :x, :z), zbrain_units)
size(h2b_zbrain)

######
## Register Zbrain to Zseries using ANTs
######

# we use histogram matching in ANTs (although maybe our method is better...?)
fixed = zseries
moving = h2b_zbrain
moving = reverse(moving,dims=2); # face right
moving = reverse(moving,dims=3); # bottom to top

# affine
# zbrain_registered = ants_register(fixed, moving)
warps = glob("*SyN_Warped.nii.gz", fishDir)
if length(warps)==0
    cmd = ants_register(fixed, moving; interpolation = "WelchWindowedSinc",
        histmatch = 0, sampling_frac = 0.25, maxiter = 200, threshold=1e-8,
        use_syn = true, synThreshold = 1e-7, synMaxIter = 200, dont_run = true)
else
    cmd = "using cached ANTs for zbrain registration"
end
cmd
## read in SyN..
zbrain_warpedname = glob("*SyN_Warped.nii.gz", fishDir)[1]
zbrain_registered = niread(zbrain_warpedname);
adj_zbrain_registered = adjust_histogram(imadjustintensity(zbrain_registered), Equalization());
## visualize results as sanity check
im = RGB.(adj_zbrain_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries
##
im[:,:,75]

## read multiMAP
oir_dir = "/data/dlab/b115/2021-07-14_rsChRmine_h2b6s_5dpf/fishfrom_2021-07-13_rschrmine_h2b6s/fish1"
# oir_920_file = joinpath(oir_dir, "multimap_zseries_920nm_ch4-gad405_ch1-sert647_chr-gcamp_take2.oir")
# oir_820_file = joinpath(oir_dir, "multimap_zseries_ch4-gad405_ch1-sert647_chr-gcamp_take2.oir")
oir_920_file = joinpath(oir_dir, "multimap_zseries_920nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom.oir")
oir_820_file = joinpath(oir_dir, "multimap_zseries_820nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom-take2.oir")
x_um, y_um, z_um = read_oir_units(oir_920_file)
@assert z_um == 1.0 # make sure not a eg 10um, quick-pan goof
# read_oir_units(oir_920_2x_file)

"Read a 4-channel image"
function read_multimap(btf_file)
    btf = collect(load(btf_file))
    oZ = size(btf, 3)
    # only blue for unknown reasons...
    btf = reshape(btf,
        size(btf)[1:2]..., 4, 2, Int(oZ//8))
    # ch1 = 646nm (red)
    channelview(btf)[1,:,:,1,:,:] = copy(channelview(btf)[3,:,:,1,:,:])
    channelview(btf)[2:3,:,:,1,:,:] .= 0
    # ch3 = 920nm (green)
    channelview(btf)[2,:,:,3,:,:] = copy(channelview(btf)[3,:,:,3,:,:])
    channelview(btf)[[1,3],:,:,3,:,:] .= 0
    #ch4 = 405nm (blue), already okay
    btf
end

# TODO why does this crash the first run only...?
# added a copy, does this fix?
@warn "need to run twice for unknown reasons to avoid error."
try
    oir_920_img = read_multimap(replace(oir_920_file, "oir" => "ome.btf"));
catch
    println("trying again...")
    oir_920_img = read_multimap(replace(oir_920_file, "oir" => "ome.btf"));
end
oir_820_img = read_multimap(replace(oir_820_file, "oir" => "ome.btf"));
##
z = 76
0.6 * Gray.(imadjustintensity(oir_920_img[:,:,3,1,z])) +
    adj_gamma_rgb(imadjustintensity(oir_820_img[:,:,4,1,z]), gamma=0.5) +
    adj_gamma_rgb(imadjustintensity(oir_820_img[:,:,1,1,z],(0,0.5)), gamma=0.5)
    # adj_gamma_rgb(imadjustintensity(oir_820_img[:,:,4,1,z]), gamma=0.2) +
    # adj_gamma_rgb(imadjustintensity(oir_820_img[:,:,1,1,z],(0,0.5)), gamma=0.2)

##
warps = glob("*SyN_Warped.nii.gz", oir_dir)
if length(warps) == 0
    moving = AxisArray(oir_920_img[:,:,3,1,:], (:y, :x, :z), (y_um, x_um, z_um));
    cmd = ants_register(fixed, moving; interpolation = "WelchWindowedSinc",
        histmatch = 0, sampling_frac = 0.25, maxiter = 200, threshold=1e-8,
        use_syn = true, synThreshold = 1e-7, synMaxIter = 200,
        save_dir=oir_dir, dont_run = true)
else
    cmd = "using cached ANTs for multiMAP registration"
end
cmd
## read registration results for 920nm multiMAP zseries
# warpedname = glob_one_file("*SyN_Warped.nii.gz", oir_dir)
warpedname = glob("*SyN_Warped.nii.gz", oir_dir)[2]
multimap_registered = niread(warpedname);
adj_multimap_registered = adjust_histogram(imadjustintensity(multimap_registered), Equalization());
## visualize results as sanity check
im = RGB.(adj_multimap_registered)
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries;
save(joinpath(plotDir, "multimap_registered_zseries.png"), im[:,:,52])
im[:,:,52]

## transform registration results for 820nm multiMAP zseries
@warn "hardcode multimap_zseries_920nm_ch4-gad405_ch1-sert647_chr-gcamp_2x-zoom for ANTs transform..."
transformPath1 = glob_one_file("20210729T154219098Z*GenericAffine.mat", oir_dir)
transformPath2 = glob_one_file("20210729T154219098Z*SyN_1Warp.nii.gz", oir_dir)
mm820_registered = zeros(size(zseries)[1:2]...,4,size(zseries,3))
# TODO: perhaps could cache this...
@threads for c in 1:4
    moving = AxisArray(oir_820_img[:,:,c,1,:], (:y, :x, :z), (y_um, x_um, z_um));
    mm820_registered[:,:,c,:] = antsApplyTransforms(fixed, moving, transformPath1, transformPath2)
end

##
adj_mm820_registered = imadjustintensity(mm820_registered);
im = RGB.(mm820_registered[:,:,3,:])
channelview(im)[2:3,:,:,:] .= 0
channelview(im)[2,:,:,:] .= adj_zseries;
##
z = 50
# im[:,:,z]
im[500:700,500:700,50]
# Gray.(adj_multimap_registered[500:700,500:700,z])
# Gray.(mm820_registered[500:700,500:700,1,z])
# Gray.(adj_zseries[500:700,500:700,z])
##
z = 20
c = 4
# mm_plane = RGB.(imadjustintensity((registered .* (registered .> 0))[:,:,z]))
mm_plane = RGB.(adjust_histogram(imadjustintensity(mm820_registered[:,:,c,z]),
    GammaCorrection(0.1)))
channelview(mm_plane)[1:2,:,:,:] .= 0

zseries_plane = Gray.(imadjustintensity(zseries[:,:,z]))
channelview(zseries_plane)[[1,3],:,:,:] .= 0
0.5*zseries_plane + mm_plane
# 0.6 * Gray.(imadjustintensity(zseries[:,:,z])) + 
#     mm_plane
##
(zseries_plane + mm_plane)[500:700,500:700]

## find imaging plane
tseries_avg = avgTseries(tseriesdir,every=10);
# TODO: this is working just okay...
# need to adjust gamma or imadjustintensity or ..?
## use xml to lookup imaging plane correspondence to Zseries

# xml = read_xml(xmlPath);
datafolders = splitpath(zseriesdir)
xmlPath = joinpath(datafolders..., datafolders[end] * ".xml")
zseries_xml = read_xml(xmlPath);
zseries_zaxes = read_all_zaxis(zseries_xml)

datafolders = splitpath(tseriesdir)
xmlPath = joinpath(datafolders..., datafolders[end] * ".xml")
tseries_xml = read_xml(xmlPath);
tseries_zaxes = read_first_zaxis(tseries_xml)

# Sequence



## This doesn't work reliably...
function find_imaging_planes(zseries, avg_imaging)
    nZ = size(zseries, 3)
    nI = size(avg_imaging,3)
    mi = zeros(nZ,nI)
    # can't use threads due to PyCall;
    # TODO: reimplement mutual_information using pure julia..?
    @showprogress for z in 1:nZ
        for i in 1:nI
            zplane = convert(Array{Float64},collect(zseries[:,:,z]))
            zplane = imresize(zplane, size(avg_imaging)[1:2])
            img = avg_imaging[:,:,i]
            # img = adjust_histogram(avg_imaging[:,:,i], GammaCorrection(0.3))
            mi[z,i] = mutual_information(zplane, img)
        end
    end
    matching_zplanes = Tuple.(argmax(mi,dims=1)[1,:])
    map(x->x[1], matching_zplanes)

end

matching_zplanes = find_imaging_planes(zseries, tseries_avg)

# ideally all three are 1:1 BUT probably spectral shift...
# sloooooow... and fails!!
# mm920_820_matching = [34, 12, 11, 11, 12, 11, 18, 18, 18, 18, 18, 18, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 28, 31, 31, 31, 31, 31, 33, 34, 34, 34, 37, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 46, 46, 50, 50, 50, 50, 50, 50, 50, 50, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 62, 62, 56, 56, 56, 50, 50, 39, 39, 39, 39, 39, 39, 39, 50, 50, 50, 39, 50, 50, 50, 39, 39, 39, 50, 39, 50, 50, 50, 50, 50, 63, 61, 63, 63, 63, 63, 63, 61, 61, 63, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 61, 62, 26, 26, 39, 40, 66, 59, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

mm920_820_matching = find_imaging_planes(mm820_registered[:,:,3,:], multimap_registered)
mm820_zseries_matching = find_imaging_planes(mm820_registered[:,:,3,:], zseries)
mm920_zseries_matching = find_imaging_planes(multimap_registered, zseries)


## read tseries avg_stim & other needed info
@warn "hardcode path for ty.h5"

expName = splitpath(tseriesdir)[end]
# avg_stim_h5_path = glob_one_file(expName*"*_avgStim.h5", fishDir)
avg_stim_h5_path = "/mnt/deissero/users/tyler//b115/2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061lstm_divide8192_avgStim_lstm.h5"
have_avg_stim_h5 = isfile(avg_stim_h5_path)
avgStim = h5read(avg_stim_h5_path, "/block1");

datafolders = splitpath(tseriesdir)
xmlPath = joinpath(datafolders..., datafolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
Z = size(etlVals,1)

slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)

target_groups = []
group_stim_freq = []
for mat in matread.(findMatGroups(slmExpDir))
    push!(target_groups, mat["cfg"]["maskS"]["targets"][1])
    push!(group_stim_freq, getMatStimFreq(mat))
end

nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)

tylerSLMDir = joinpath(fishDir, "slm")
slmTxtFile = regex_glob(r".*(?<!trialOrder)\.txt$", tylerSLMDir)
@assert length(slmTxtFile) == 1 slmTxtFile # if not, need to be careful to choose
slmTxtFile = slmTxtFile[2]
##
stimGroupDF = CSV.File(open(read, slmTxtFile), header=["filepath", "powerFraction"]) |> DataFrame
stimGroupDF = stimGroupDF[trialOrder,:]

mat = matread.(findMatGroups(slmExpDir)[1])
max_targets_center = maximum(map(mat->maximum(maximum.(matread(mat)["cfg"]["exp"]["maskS"]["targetsCenter"])), findMatGroups(slmExpDir)))
@assert  max_targets_center < 512 "bad targets center value, don't trust results..."
slmNum = getSLMnum(mat)


is1024 = size(avgStim,1)==1024
tseriesH, tseriesW = size(avgStim)[1:2]
tseries_lateral_unit = read_microns_per_pixel(xmlPath)[1]

zOffset = getzoffset(expDate, slmNum)
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
    is1024=is1024, zOffset=zOffset)


voltageFiles = glob("*VoltageRecording*.csv", tseriesdir)
stimStartIdx, stimEndIdx, frameStartIdx = getStimTimesFromVoltages(voltageFiles[1], Z);
cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
cells[!, :stimFreq] = map(g->group_stim_freq[trialOrder][g], cells.stimGroup)
# cells[!, :laserPower] = round.(typeof(1.0mW),
#     map(g->stimGroupDF.powerFraction[g], cells.stimNum) .* slmpowerPerCell, digits=1)

targetSizePx = spiral_size(expDate, tseries_lateral_unit)

## make figure
@warn "hardcode volRate, max_frames, nStimuli, Z"
volRate = 3
max_frames = 15
nStimuli = 3
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)

imap_cmap = transparent_cmap(plt.cm.Reds, max_alpha=0.7)
bg_cmap = transparent_cmap(plt.cm.gray, max_alpha=0.3)
mm4_cmap = transparent_cmap(plt.cm.Blues, max_alpha=0.5)
mm1_cmap = transparent_cmap(plt.cm.Greens, max_alpha=0.5)
# TODO: plot I-map & multiMAP on same plot

max_frames = 30
@warn "TODO: double check these"
window = minimum([Int(ceil(3*volRate)), max_frames])
# @assert (window <= post) & (window <= pre)
# @assert Z == 1
tseries_cmax = 2.5
tseries_cmin = -0.5

# tseries_cmax = 1.0
tseries_cmax = 0.8
tseries_cmin = 0.0

mm4_cmax = quantile(mm820_registered[:,:,4,:][:], .999)
mm4_cmin = 0.0
mm1_cmax = quantile(mm820_registered[:,:,1,:][:], .999)
mm1_cmin = 0.0

# tseries_cmax = 0.5
# tseries_cmin = -0.1
# tseries_cmax = 4
# tseries_cmin = -0.75
# tseries_cmin = -4
# cnorm = matplotlib.colors.TwoSlopeNorm(vmin=tseries_cmin,vcenter=0,vmax=cmax)

recording_folder = splitpath(tseriesdir)[end-2]
fish_name = splitpath(tseriesdir)[end-1]
analysis_name = "lstm-multiMAP-zoom2x"
# analysis_name = "lstm"
plotDir = joinpath(fishDir, "plots-denoised")

@warn "df_f denominator epsilon may have changed"
for stimNum in 1:nStimuli
# for stimNum in [2]
    # f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    # TODO: should window be immediately after stim..? right now it's
    # ~2 seconds after stim...?
    f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    f0 = mean(avgStim[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    df = f - f0
    # originally no epsilon prior to July 2021
    df_f = df./(f0 .+ 0.05)
    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    if Z > 5
        figW,figH = (figB*Z/2, figB*2)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(2,Int(Z/2), figsize=(figW,figH))
    else
        figW,figH = (figB*Z, figB)
        global fig = plt.figure(figsize=(figW,figH))
        fig, axs = plt.subplots(1,Z, figsize=(figW,figH))
    end
    if Z==1
        axs = [axs]
    end
    # ax = axs[2]
    # ax = axs
    for z in 1:Z
        ax = axs[z]
        # global cim = ax.imshow(df_f[:,:,z], cmap="RdBu_r",
        #     norm=cnorm)
        match_z = matching_zplanes[z]
        # fix scaling issue
        df_size = size(df_f)[1:2]
        # zseries
        ax.imshow(imresize(zseries[:,:,match_z], df_size), cmap=bg_cmap)
        # influence map
        global cim = ax.imshow(df_f[:,:,z], cmap=imap_cmap,
            clim=(tseries_cmin,tseries_cmax))
        # multiMAP
        ax.imshow(imresize(mm820_registered[:,:,4,match_z], df_size), cmap=mm4_cmap,
            clim=(mm4_cmin,mm4_cmax))
        ax.imshow(imresize(mm820_registered[:,:,1,match_z], df_size), cmap=mm1_cmap,
            clim=(mm1_cmin,mm1_cmax))
        ax.set_axis_off()
        # TODO: fix error on next line
        ax.set_title("$(Int(round(etlVals[z],digits=0)))μm")
        # this will make extra circles (1 extra per repetition...)
        for (x,y,targetZ) in eachrow(unique(cells[cells.stimNum .== stimNum,[:x,:y,:z]]))
            if z == targetZ
                circle = matplotlib.patches.Circle((x,y), targetSizePx, color="k",
                    fill=false, lw=0.4, alpha=0.3)
                ax.add_patch(circle)
            end
        end
        # previous (but not this stim) targets
        # for (x,y,z) in eachrow(unique(cells[cells.stimNum .!= stimNum,[:x,:y,:z]]))
        #     @assert z == 1
        #     circle = matplotlib.patches.Circle((x,y), targetSizePx, color="k",
        #         fill=false, lw=0.5, alpha=0.5)
        #     ax.add_patch(circle)
        # end
    end

    # axs[1].imshow(stim_masks[:, :,z,stimNum], cmap="gray")

    # cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
    # plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
    #     norm=cnorm)
    # may need to adjust if colorbar is cutoff
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.15, 0.0075, 0.7])
    # cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
    cbar = fig.colorbar(cim, cax=cbar_ax)
    path = joinpath(plotDir,"$(recording_folder)_$(fish_name)_$(expName)_$(analysis_name)_stim$stimNum")
    @show path*".svg"
    fig.savefig(path*".svg", dpi=600)
    fig.savefig(path*".png", dpi=600)
end
