# need X11 forwarding for ssh, `ssh -X osprey`, then `echo $DISPLAY`, and fill below
# ENV["DISPLAY"] = "localhost:11"
# ENV["DISPLAY"] = ":0"
using FileIO, NRRD, HDF5, MAT, Images,
    Unitful, AxisArrays, StaticArrays, CoordinateTransformations,
    TestImages, NRRD, LinearAlgebra, ImageMagick,
    LibExpat, Glob, ANTsRegistration, ImageSegmentation, PyCall,
    NIfTI, ImageContrastAdjustment, Images, #ImageView,
    ImageSegmentation, ImageMagick, Random, StatsBase,
    ImageDraw, ImageFiltering, Glob, HDF5, Arrow,
    Dates, Distributed, SharedArrays, CSV, DataFrames, Statistics, Lensman
using Gadfly
import Plots
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib
# matplotlib.use("Agg")

font = Dict("family" => "Arial",
        "weight" => "normal",
        "size"   => 10)

matplotlib.rc("font", family="Arial", weight="normal", size=10)

SMALL_SIZE = 10
matplotlib.rc("font", size=SMALL_SIZE)          # controls default text sizes
matplotlib.rc("axes", titlesize=SMALL_SIZE)     # fontsize of the axes title
matplotlib.rc("axes", labelsize=SMALL_SIZE)    # fontsize of the x and y labels
matplotlib.rc("xtick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc("ytick", labelsize=SMALL_SIZE)    # fontsize of the tick labels
matplotlib.rc("legend", fontsize=SMALL_SIZE)    # legend fontsize
matplotlib.rc("figure", titlesize=SMALL_SIZE)  # fontsize of the figure title

matplotlib.rcParams["figure.dpi"] = 600
matplotlib.rcParams["savefig.dpi"] = 600

using Base.Iterators: peel
import Unitful: μm
##
tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
# tseriesRootDir = "/mnt/deissero/users/tyler/b115"
# tseriesRootDir = "/data/dlab/b115"



# tseriesRootDir = "/mnt/deissero/users/tyler/b115"
# tseriesRootDir = "/scratch/b115"
# tifDir = "$tseriesRootDir/2021-02-16_6f_h33r_f0_6dpf/fish2/TSeries-lrhab-raphe-control-129trial-052/"
# tifDir = "$tseriesRootDir/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-039/"
# tifDir = "$tseriesRootDir/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish2/TSeries-lrhab_raphe_40trial-045/"
# tifDir = "$tseriesRootDir/2020-11-02_elavl3-chrmine-Kv2.1_h2b6s_5dpf/fish1/TSeries-lrhab_raphe_40trial-040/"

# tifDir = "$tseriesRootDir/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1/TSeries-lrhab_raphe_40trial-023"
tifDir = "$tseriesRootDir/2020-10-26_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish2/TSeries-lrhab_raphe_stim-40trial-proper-zoffset-034"
Z = 5

# tifDir = "$tseriesRootDir/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-lrhab_raphe_stim-40trial-038"
# tifDir = "$tseriesRootDir/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries-lrhab_raphe_stim-40trial-038"

# slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
# slmDir = "$tseriesRootDir/tyler/b115/SLM_files/"
# slmDir = "/oak/stanford/groups/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment"
slmDir = "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files"
# slmDir = "$tseriesRootDir/SLM_files/"

expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
recording_folder = splitpath(tifDir)[end-2]
fish_name = splitpath(tifDir)[end-1]

plotDir = "/oak/stanford/groups/deissero/users/tyler/plots/2021_chrmine-structure"
# plotDir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure"
# plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end
plotDir
##
tseries = loadTseries(tifDir);
(H, W, Z, T) = size(tseries)
##
voltageFile = joinpath(tifDir, expName*"_Cycle00001_VoltageRecording_001.csv")
if ~isfile(voltageFile)
    # avoid scanning list of many files...
    voltageFile = glob("*VoltageRecording*.csv", tifDir)[1]
end
stimStartIdx, stimEndIdx = getStimTimesFromVoltages(voltageFile, Z)
# @assert length(stimStartIdx) == 129
@assert length(stimStartIdx) == 120 "found only $(length(stimStartIdx)) stimuli"

## single stim example
# imshow(imadjustintensity(tseries[:,:,:,stimStartIdx[1]-1:stimEndIdx[1]+1]))

## read slm stim files
dataFolders = splitpath(tifDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
volRate = frameRate / Z

pre = Int(ceil(5*volRate))
post = Int(ceil(5*volRate))

# if imaging many planes, may need to read another xml file since didn't get all planes in first file
# etlVals = etlVals[1:end-1]
@assert length(etlVals) == Z
@info "assume imaging from top-down"

#
slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)


nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
@show nTrials, nStimuli
@assert nTrials == size(stimStartIdx,1) # check TTL start-times match 

target_groups = [mat["cfg"]["maskS"]["targets"][1]
    for mat in matread.(findMatGroups(slmExpDir))]

# Array of matrices indicating x, y z (or was it y, x, z..?)
groupLocs = mapTargetGroupsToPlane(target_groups, etlVals, is1024=false)
groupLocs = map(x->Int.(round.(x)), groupLocs)

# stimDuration = Int(ceil(2*volRate))

##
microscope_units = (0.6299544139175637μm, 0.6299544139175637μm, 2.0μm)
if W == 1024
    microscope_units = microscope_units ./ 2
end

targetSizePx = (7μm * 14.4/25) / microscope_units[1]
regionMasks = constructGroupMasks(groupLocs, H, W, Z, targetSizePx=targetSizePx);

@assert all(sum(regionMasks,dims=[1,2,3]) .> 0)

## acquire in-region stim traces
winSize = Int(ceil(volRate*2))
before = Int(ceil(volRate*5))
after = Int(ceil(volRate*14))
region_trace_df = DataFrame(region=Int[], time=Float64[], df_f=Float64[], stimIdx=Int64[])
for (i, stim_start, stim_stop, stim) in zip(1:length(trialOrder), stimStartIdx, stimEndIdx, trialOrder)
    s = stim_start-before
    e = stim_stop+after
    time = collect(s:e) .- stim_start * 1.0
    time ./= volRate
    region_trace = extractTrace(tseries[:,:,:,s:e], regionMasks[:,:,:,stim])
    region_trace = imageJkalmanFilter(region_trace)
    
    f0 = mean(region_trace[1:winSize])
    df_f = @. (region_trace - f0) / (f0 + 10)
    region_trace_df = vcat(region_trace_df,
        DataFrame(region=stim, time=time, df_f = df_f, stimIdx=i))
end
region_trace_df
open(joinpath(fishDir, "$(expName)_fluorescence.arrow"), "w") do io
    Arrow.write(io, region_trace_df)
end

"wrote fluorescence arrow!"
## plot the traces!

chrmine_paper = Gadfly.Theme(
    line_width=1mm,
    minor_label_font="Arial",
    major_label_font="Arial",
    key_title_font="Arial",
    key_label_font="Arial",
    major_label_font_size=10pt,
    minor_label_font_size=10pt,
    key_title_font_size=10pt,
    key_label_font_size=10pt
)

meanDF = combine(groupby(region_trace_df, ["time", "region"]),
    :df_f => mean, :df_f => maximum, :df_f => minimum)
    # :df_f => mean, :df_f => max_95ci, :df_f => min_95ci)

##

line = plot(meanDF, x=:time,
    xgroup=:region,
    y="df_f_mean",
    # ymin="df_f_minimum",
    # ymax="df_f_maximum",
    # ymin="df_f_min_95ci",
    # ymax="df_f_max_95ci",
    chrmine_paper,
    # Geom.subplot_grid(Geom.line, Geom.ribbon))
    Geom.subplot_grid(Geom.line))
path = joinpath(plotDir,"$(expName)_region_trace.svg")
Gadfly.draw(SVG(path, 12cm, 6cm), line)
@show path
line
## avg stim effect

# warmup JIT
# @time trialAverage(tseries, stimStartIdx[1:3], stimEndIdx[1:3], trialOrder[1:3],
#     pre=pre, post=post);
avgStim = trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder;
    pre=pre, post=post);
##
try
    h5write(joinpath(fishDir,expName*"_avgStim.h5"), "/block1", avgStim)
catch
end
##
# avgStim = h5read(joinpath(fishDir,expName*"_avgStim.h5"), "/block1");
(H, W, Z, nStim) = size(avgStim)
## STA

figB = 1.6
figW,figH = (figB*5, figB)

window = Int(ceil(3*volRate))
@assert (window < post) & (window < pre)
# cmax = 2.5
# cmin = -0.5
cmax = 4
cmin = -0.75
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
# cnorm = matplotlib.colors.DivergingNorm(vmin=cmin,vcenter=0,vmax=cmax)

for stimNum in 1:nStimuli
    f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    f0 = mean(avgStim[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    df = f - f0
    global df_f = df./f0

    # cmax = percentile(df_f[:],99.9)
    # cmin = percentile(df_f[:],0.1)
    # global fig = plt.figure(figsize=(figW,figH))
    global fig, axs = plt.subplots(1,Z, figsize=(figW,figH))
    for (z,ax) in enumerate(axs)
        # plt.axis("off")
        global cim = ax.imshow(df_f[:,:,z], cmap="RdBu_r",
            norm=cnorm)
        ax.set_axis_off()
        ax.set_title("$(Int(round(etlVals[z],digits=0)))μm")
    end
    # cmax = percentile(abs.(df_f[:,:,1][:]),99.9)
    # plt.imshow(hcat([df_f[:,:,z] for z in 1:Z]...), cmap="RdBu_r",
    #     norm=cnorm)
    fig.subplots_adjust(right=0.96)
    cbar_ax = fig.add_axes([0.97, 0.15, 0.0075, 0.7])
    # cbar = fig.colorbar(cim, ticks=[0,1,2], cax=cbar_ax)
    cbar = fig.colorbar(cim, cax=cbar_ax)
    path = joinpath(plotDir,"$(recording_folder)_$(fish_name)_$(expName)__stim$stimNum.svg")
    @show path
    fig.savefig(path, dpi=600)
end
fig
