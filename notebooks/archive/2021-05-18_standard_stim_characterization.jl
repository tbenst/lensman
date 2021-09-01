## compare 
ENV["DISPLAY"] = "localhost:11.0"
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
using ImageView
import Gadfly
using Unitful: μm, m, s, mW
import Base.Threads.@threads
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib

##
ON_SHERLOCK = read(`hostname`,String)[1:2] == "sh"
if ON_SHERLOCK
    tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
else
    # tseriesRootDir = "/data/dlab/b115"
    tseriesRootDir = "/scratch/b115"
    # tseriesRootDir = "/mnt/deissero/users/tyler/b115"
end

# newer
if ON_SHERLOCK
    slmDir = "/oak/stanford/groups/deissero/users/tyler/b115/SLM_files"
else
    slmDir = "/mnt/deissero/users/tyler/b115/SLM_files"
    # slmDir = "/mnt/deissero/users/tyler/slm/mSLM/SetupFiles/Experiment"
    # slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
    # older
    # slmDir = "/mnt/deissero/users/tyler/b115/SLM_files/"
    # slmDir = "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/"
end


# tseriesDir = "/data/dlab/b115/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-1024cell-32concurrent-4freq-054"
# tseriesDir = "/data/dlab/b115/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-256cell-8concurrent-4freq-055"
# can't fit in memory on lensman, so use deis
# tseriesDir = joinpath(tseriesRootDir, "2021-02-16_6f_h33r_f0_6dpf/fish2/TSeries-256cell-8concurrent-4freq-051")
# tseriesDir = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043")
# tseriesDir = joinpath(tseriesRootDir, "2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2/TSeries-128cell-4concurrent-3power-skip7-044")
# tseriesDir = "$tseriesRootDir/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043"
# tseriesDir = "$tseriesRootDir/2021-01-26_rsChRmine_6f_7dpf/fish2/TSeries-32concurrent-256trial-2rep-4power-045" # only one stim...?
# tseriesDir = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-10zplane-077"
# tseriesDir = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-076"
# tseriesDir = "$tseriesRootDir/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060"

# tseriesDir = "$tseriesRootDir/2021-03-16_h33r-chrmine_h2b6s/fish4/TSeries_64cell_8concurrent_2power_8rep-607"
# tseriesDir = "$tseriesRootDir/2021-03-16_rschrmine_h2b6s/fish3/TSeries_64cell_8concurrent_2power_8rep-407"
# tseriesDir = "$tseriesRootDir/2021-03-16_h2b6s/fish1/TSeries_64cell_8concurrent_2power_8rep-207"
# tseriesDir = "$tseriesRootDir/2021-03-16_wt-chrmine_h2b6s/fish2/TSeries_64cell_8concurrent_2power_8rep-221"

# tseriesDir = joinpath(tseriesRootDir, "2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2/TSeries-128cell-4concurrent-3power-skip7-044")
# tseriesDir = "$tseriesRootDir/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043"
# tseriesDir = "$tseriesRootDir/2021-01-26_rsChRmine_6f_7dpf/fish2/TSeries-32concurrent-256trial-2rep-4power-045" # only one stim...?
# tseriesDir = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-10zplane-077"
# tseriesDir = "$tseriesRootDir/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060"
# tseriesDir = "$tseriesRootDir/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060"
# tseriesDir = "$tseriesRootDir/2021-01-26_rsChRmine_6f_7dpf/fish2/TSeries-1024cell-32concurrent-048"
# tseriesDir = "$tseriesRootDir/2021-03-09_wt-chrmine-gc6f/fish1/TSeries-32cell-8concurrent-10MHz-8rep-065"
# tseriesDir = "$tseriesRootDir/2021-03-09_wt-chrmine-gc6f/fish1/TSeries-32cell-8concurrent-10MHz-8rep-065"


# compare to old 1024 stim experiment with 90% of cells responding (but some plateaus...)
# tseriesDir = "$tseriesRootDir/2020-12-15_h2b6s_chrmine_kv2.1_5dpf/fish2/TSeries-1024cell-32concurrent-038"
# tseriesDir = "$tseriesRootDir/2020-12-17_h2b6s_chrmine-kv2.1_7dpf/fish1/TSeries_n64_b2_r8-044"
# tseriesDir = "$tseriesRootDir/2020-11-30_elavl3-chrmine-Kv2.1_h2b6s_6dpf/fish1/fov1-tectum/TSeries-64groups-4cell-stim-069"

# tseriesDir = "$tseriesRootDir/2021-03-30_wt-chrmine_6dpf_h2b6s/fish1/TSeries-28cell-1concurrent-2power-10trial-052"
# tseriesDir = "$tseriesRootDir/2021-04-13_wt-chrmine_6dpf_h2b6s/fish1/TSeries-16cell-50rep-actually-stim-147"
# tseriesDir = "$tseriesRootDir/2021-04-19_wt-chrmine_5dpf_6f/fish1/TSeries-15cell-5concurrent-5rep-actual-trigger-005"
# tseriesDir = "$tseriesRootDir/2021-04-19_wt-chrmine_5dpf_6f/fish1/TSeries-15cell-5concurrent-5rep-actual-trigger-005"
# tseriesDir = "$tseriesRootDir/2021-04-20_h33r-chrmine_6dpf_6f/fish2/TSeries-20cell-1concurrent-20trial-047"

# tseriesDir = "$tseriesRootDir/2021-05-18_hsChRmine_h2b6s_6dpf/fish4/TSeries-lrhab-control-91trial-4Mhz-044"
# tseriesDir = "$tseriesRootDir/2021-05-18_rsChRmine_h2b6s_6dpf/fish1/TSeries-lrhab-control-91trial-4Mhz-047"
# tseriesDir = "$tseriesRootDir/2021-05-18_rsChRmine_h2b6s_6dpf/fish5/TSeries-lrhab-control-91trial-4Mhz-045"
# tseriesDir = "$tseriesRootDir/2021-05-18_rsChRmine_h2b6s_6dpf/fish5/TSeries-35cell-20rep-40s-dark-4Mhz-059"

# tseriesDir = "$tseriesroot/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-lrhab-118trial-061"
# tseriesDir = "$tseriesRootDir/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
# tseriesDir = "$tseriesRootDir/2021-06-02_rsChRmine-h2b6s/fish2/TSeries-titration-192trial-062"
# tseriesDir = "$tseriesRootDir/2021-06-08_rsChRmine-h2b6s/fish2/TSeries-titration-192trial-062"
# tseriesDir = "$tseriesRootDir/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-1232021-06-21_6pm"
tseriesDir = "$tseriesRootDir/2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
# debug by looking at 3region..?
# tseriesDir = "$tseriesRootDir/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish2/TSeries_lrhab_raphe_40trial-044/"

##
# tseriesDir = "$tseriesRootDir/2021-03-16"
# gcampDir = "/scratch/b115/2021-03-16_h2b6s/fish1"
# h33rDir = "/scratch/b115/2021-03-16_h33r-chrmine_h2b6s/fish4"
# rsDir = "/scratch/b115/2021-03-16_rschrmine_h2b6s/fish3/"
# wtDir = "/scratch/b115/2021-03-16_wt-chrmine_h2b6s/fish2/"

# tseriesDir = "$tseriesRootDir/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060"
# possibly compare to...
# 2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/ (4power)
# 2021-02-15_wt_chrmine_gc6f/fish1/TSeries-1024cell-4freq-skip-first-066 (4freq; too large for memory on lensman)


# analysis_name = "lstm_denoise_only"
analysis_name = ""
# analysis_name = "kalman"

tyh5Path = tseriesDir * ".ty.h5"

if occursin("freq", tseriesDir)
    exp_param = :stimFreq
elseif occursin("power", tseriesDir)
    exp_param = :laserPower
end

if tseriesDir[end] == "/"
    tseriesDir = tseriesDir[1:end-1]
end


if occursin("freq", tseriesDir)
    exp_param = :stimFreq
elseif occursin("power", tseriesDir)
    exp_param = :laserPower
end

if tseriesDir[end] == "/"
    tseriesDir = tseriesDir[1:end-1]
end

fishDir = joinpath(splitpath(tseriesDir)[1:end-1]...)
expName = splitpath(tseriesDir)[end]
recording_folder = splitpath(tseriesDir)[end-2]
fish_name = splitpath(tseriesDir)[end-1]
tylerSLMDir = joinpath(fishDir, "slm")

tyh5Path = joinpath(fishDir,"TSeries-lrhab-titration-1232021-06-21_6pm" * ".ty.h5")
# tyh5Path = tseriesDir*".ty.h5"
##
"Take Float32 [0,1], and convert to UInt16"
function float2uint(x)
    # PV max is 8192..? or at least, that's what we used in babelfish
    ret = round.(x .* 8192, digits=0)
    clamp!(ret,0,2^16-1)
    convert(Array{UInt16},ret)
end

# TODO: lazy read hdf5 wrapper...?

function read_tyh5(tyh5_path)
    @info "assume WHZCT, convert to HWZT"
    h5 = h5open(tyh5Path,"r")
    # dset = h5["/imaging/per_pixel_lstm_denoised"]
    # dset = h5["/imaging/per_pixel_lstm_denoised_maybe_longer_time"]
    dset = h5["/imaging/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm"]
    println(ndims(dset))
    if ndims(dset) == 5
        @assert size(dset,4)==1
        W,H,Z,C,T = size(dset)
        dtype = eltype(dset)
        if dtype == Float32
            # convert to UInt16, asume max value of 1...
            f_convert = float2uint
        else
            f_convert = x -> x
        end
        tseries = zeros(UInt16, H, W, Z, T)
        # chunk by z for memory efficiency..
        @showprogress for z=1:Z
            # drop singleton channel
            dat = f_convert(dset[:,:,z,1,:])
            dat = permutedims(dat, (2,1,3))
            tseries[:,:,z,:] = dat
        end
    else
        error("stop")
        tseries = h5read(tyh5Path, "/imaging/per_pixel_lstm_denoised")
        # tseries = h5read(tyh5Path, "/imaging/raw")
        # @show size(tseries)
        # @assert size(tseries,2)==1
        # tseries = tseries[:,1,:,:,:]
        tseries = permutedims(tseries, (2,1,3,4))
        tseriesDir = joinpath(fishDir, expName)

        # tseries = h5read(tyh5Path, "/imaging/raw")
        # drop singleton channel
        # @assert size(tseries,4)==1
        # tseries = permutedims(tseries, (2,1,3,4,5))
        # tseries = tseries[:,:,:,1,:];
        # tseriesDir = joinpath(fishDir, expName)
    end
    tseries
end

if isfile(tyh5Path)
    @info "using ty.h5 file"
    tseries = read_tyh5(tyh5Path)
else
    tseries = loadTseries(tseriesDir);
end;
##

GC.gc()
##
(H, W, Z, T) = size(tseries)
if Z > 1000
    tseries = reshape(tseries, H, W, 1, Z*T);
    (H, W, Z, T) = size(tseries)
end
@show (H, W, Z, T)
##
plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end

# TODO: wrap more of cell in function to reduce inline code..?
# for 2021-03-16, can read from oak...
# voltageFile = glob("*VoltageRecording*.csv", tseriesDir)[1]
voltageFiles = glob("*VoltageRecording*.csv", tseriesDir)

# read slm stim files
dataFolders = splitpath(tseriesDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
Z = size(etlVals,1)
volRate = frameRate / Z
try
    global trialOrder, slmExpDir
    slmExpDir = joinpath(slmDir,Dates.format(expDate, "dd-u-Y"))
    trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)
catch
    # not sure if needed
    # hack in case matlab code uses date program started not experiment run
    slmExpDir = joinpath(slmDir,Dates.format(expDate - Day(1), "dd-u-Y"))
    @info "try previous day"
    trialOrder, slmExpDir = getTrialOrder(slmExpDir, expDate)
end

## read power
try
    global slmTxtFile
    slmTxtFile = regex_glob(r".*(?<!trialOrder)\.txt$", tylerSLMDir)
catch
    global slmTxtFile
    slmTxtFile = regex_glob(r".*(?<!trialOrder)\.txt$", fishDir)
end
@assert length(slmTxtFile) == 1 slmTxtFile # if not, need to be careful to choose
slmTxtFile = slmTxtFile[1]

##
stimGroupDF = CSV.File(open(read, slmTxtFile), header=["filepath", "powerFraction"]) |> DataFrame
stimGroupDF = stimGroupDF[trialOrder,:]

mat = matread.(findMatGroups(slmExpDir)[1])
max_targets_center = maximum(map(mat->maximum(maximum.(matread(mat)["cfg"]["exp"]["maskS"]["targetsCenter"])), findMatGroups(slmExpDir)))
@assert  max_targets_center < 512 "bad targets center value, don't trust results..."
slmNum = getSLMnum(mat)
zOffset = getzoffset(expDate, slmNum)
##
firstTargetGroup = matread.(findMatGroups(slmExpDir)[1])
powerPerCell = firstTargetGroup["cfg"]["mode"]["BHV001"]["FOV"]["PowerPerCell"]
slm1Power, slm2Power = slmpower(expDate)
if slmNum == 1
    slmpowerPerCell = slm1Power * powerPerCell / 1000
elseif slmNum == 2
    slmpowerPerCell = slm2Power * powerPerCell / 1000
end

if length(voltageFiles) == 1
    stimStartIdx, stimEndIdx, frameStartIdx = getStimTimesFromVoltages(voltageFiles[1], Z);
else
    @warn "HARDCODED TperTrial"
    TperTrial = Int(T / 20)
    stimStartIdx, stimEndIdx, frameStartIdx = getStimTimesFromVoltages(voltageFiles, Z, TperTrial);
end

if length(voltageFiles) == 1
    allWithin(diff(stimStartIdx),0.05)
    # adjust for expected number of stimuli....
    # @assert length(stimStartIdx) == 32
    @show length(stimStartIdx)
    allWithin(stimEndIdx .- stimStartIdx,0.05)
    avgStimDuration = mean(stimEndIdx .- stimStartIdx)
    # adjust for volume rate / stim duration...
    @assert 30 > avgStimDuration > 3
else
    @warn "multiple voltage files so not checking asserts in case a trial structure..?"
end

# Assumes no sequence of stim
# for 5Hz clock

target_groups = []
group_stim_freq = []
for mat in matread.(findMatGroups(slmExpDir))
    push!(target_groups, mat["cfg"]["maskS"]["targets"][1])
    push!(group_stim_freq, getMatStimFreq(mat))
end
##
nStimuli = maximum(trialOrder)
nTrials = size(trialOrder,1)
# nTrialsPerStimulus = Int(size(trialOrder,1) / nStimuli)
# @assert nTrialsPerStimulus == 1
# cells = Array{Array{Float64,1}}
is1024 = size(tseries,1)==1024
targetsWithPlaneIndex = mapTargetGroupsToPlane(target_groups, etlVals,
is1024=is1024, zOffset=zOffset)
targetsWithPlaneIndex = map(x->Int.(round.(x, digits=0)), targetsWithPlaneIndex)

cells = makeCellsDF(targetsWithPlaneIndex, stimStartIdx, stimEndIdx, trialOrder)
cells[!, :stimFreq] = map(g->group_stim_freq[trialOrder][g], cells.stimGroup)
cells[!, :laserPower] = round.(typeof(1.0mW),
    map(g->stimGroupDF.powerFraction[g], cells.stimNum) .* slmpowerPerCell, digits=1)

lateral_unit = microscope_lateral_unit(W)
targetSizePx = spiral_size(expDate, lateral_unit)

## influence maps
max_frames = stimStartIdx[1]-1 
max_time = (max_frames-1)/volRate # otherwise bad index issues
if volRate > 10
    nseconds = 3
else
    nseconds = 5
    nseconds = minimum([nseconds, max_time])
end
pre = Int(ceil(nseconds*volRate))+1
post = Int(ceil(nseconds*volRate))+1

# avg_stim_h5_path = joinpath(fishDir,expName*"$(analysis_name)_avgStim.h5")
avg_stim_h5_path = joinpath(fishDir,expName*"$(analysis_name)_avgStim_lstm.h5")
have_avg_stim_h5 = isfile(avg_stim_h5_path)
# have_avg_stim_h5 = false # temp force refresh
if have_avg_stim_h5
    avgStim = h5read(avg_stim_h5_path, "/block1");
else
    avgStim = trialAverage(tseries, stimStartIdx, stimEndIdx, trialOrder;
        pre=pre, post=post);
    h5write(avg_stim_h5_path, "/block1", avgStim)
end;

##
figB = 1.6
# 2 for extra stim mask
# figW,figH = (figB*1.1, figB)


window = minimum([Int(ceil(3*volRate)), max_frames])
@assert (window <= post) & (window <= pre)
# @assert Z == 1
cmax = 2.5
cmin = -0.5
# cmax = 4
# cmin = -0.75
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)

for stimNum in 1:nStimuli
    f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
    f0 = mean(avgStim[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
    df = f - f0
    df_f = df./f0
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
        global cim = ax.imshow(df_f[:,:,z], cmap="RdBu_r",
            norm=cnorm)
        ax.set_axis_off()
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

## kalman denoise
@time denoised = imageJkalmanFilter(tseries);
##
avg_stim_h5_path = joinpath(fishDir,expName*"$(analysis_name)_avgStim.h5")
have_avg_stim_h5 = isfile(avg_stim_h5_path)
if have_avg_stim_h5
    avgStim = h5read(avg_stim_h5_path, "/block1");
else
    avgStim = trialAverage(denoised, stimStartIdx, stimEndIdx, trialOrder;
        pre=pre, post=post);
    h5write(avg_stim_h5_path, "/block1", avgStim)
end;

##
stimLocs = map(CartesianIndex ∘ Tuple, eachrow(cells[:,[2,1]]))

stimPlane = Int(mean(map(x->mean(x[:,3]), targetsWithPlaneIndex)))
avgImage = dropdims(mean(tseries[:,:,stimPlane,:], dims=3), dims=3)
avgImageAdj = adjust_gamma(imadjustintensity(avgImage), 0.5)
avgImageAdj = RGB.(avgImageAdj)
channelview(avgImageAdj)[[1,3],:,:,:] .= 0

avgImgWithTargets = addTargetsToImage(copy(avgImageAdj), cartIdx2Array(stimLocs),
    targetSize=targetSizePx)
save(joinpath(plotDir,"$(expName)_avgImgWithTargets.png"), avgImgWithTargets)

# winSize = Int(ceil(5*volRate))
winSize = Int(ceil(2*volRate))
# delay=Int(ceil(volRate*2))
delay=0

roiMask = constructROImasks(cells, H, W, Z, targetSizePx=targetSizePx);

cellsDF = makeCellsDF(tseries, cells, roiMask, winSize=winSize, delay=delay);
# cellsDF2 = makeCellsDF(tseries, cells2, roiMask, winSize=winSize, delay=delay);

rankings = sortperm(cellsDF.df_f, lt=(>))


meanDf_f = combine(groupby(cellsDF, :cellID), :df_f => mean).df_f_mean
cellIDrankings = sortperm(meanDf_f, lt=(>))

# rankings2 = sortperm(cellsDF2.df_f, lt=(>))
numNaN = sum(isnan.(cellsDF.df_f))
@assert numNaN == 0

@show sum(cellsDF.df_f.>0.1), sum(cellsDF.df_f.>0.4), sum(cellsDF.df_f.>1.)
# @show sum(cellsDF2.df_f.>0.1), sum(cellsDF2.df_f.>0.4), sum(cellsDF2.df_f.>1.)

# only 26 are >40% for two stim in control...
# intersect(cellsDF.x[cellsDF.df_f.>0.4], cellsDF2.x[cellsDF2.df_f.>0.4])
# @show first(cellsDF[rankings[numNaN+1:end],:],50)

cellsdf_arrow_path = tseriesDir*"_cellsDF.arrow"
if ~isfile(cellsdf_arrow_path)
    open(cellsdf_arrow_path, "w") do io
        tempDF = copy(cellsDF)
        tempDF[!, :laserPower] ./= 1mW
        Arrow.write(io, tempDF)
    end
end

## sanity check that stim mask looks right
maxT = minimum([5000, size(tseries,4)])
avg_image = mean(tseries[:,:,:,1:maxT], dims=4)[:,:,1,1];
avg_image = imadjustintensity(adjust_gamma(avg_image,0.2))


stim_masks = constructGroupMasks(targetsWithPlaneIndex, H, W, Z)
Gray.(imadjustintensity(sum(stim_masks, dims=4)[:,:,1,1]))
imshow(stim_masks)

avgImgWithTargets
##
# histogram(cellsDF.df_f)
p = histogram(combine(groupby(cellsDF, [:cellID, :laserPower]), :df_f => mean).df_f_mean)
savefig(p, joinpath(plotDir,"$(expName)_df_histogram.svg"))
p
# groupbycellsDF

##
# volRate = 30
before = Int(ceil(volRate*10))
after = Int(ceil(volRate*20))
# best
plots = []
nplots = 64
##
# nplots = 8
# for idx in rankings[end:-1:end-nplots+1]
# for idx in 1:nplots
# plots = Array{Int}(undef, length(cellIDrankings))
# can we @threads Array{Int}(undef, length(cellIDrankings))
@showprogress for cellID in cellIDrankings[1:nplots]
# @threads for (i,cellID) in enumerate(cellIDrankings)
    indices = findall(cellsDF.cellID .== cellID)
    push!(plots, plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after, colorBy=exp_param))
    # plots[i] = plotStim(tseries,roiMask,cells,indices, volRate, before=before, after=after, colorBy=exp_param)
end
##
p = plot(plots..., layout=(8,8), legend=true, size=(1024*4,1024*4))
# savefig(p, joinpath(plotDir,"$(expName)_top40_traces_rollmedian.png"))
savefig(p, joinpath(plotDir,"$(expName)_top64_traces.svg"))
p
##
max_laser = maximum(cellsDF.laserPower)
cells_meanDF = combine(groupby(cellsDF[cellsDF.laserPower .== max_laser,:],
    [:cellID, :x, :y]), :df_f => mean)
cmax = maximum(abs.(cells_meanDF.df_f_mean[(~).(isnan.(cells_meanDF.df_f_mean))]))
p_df_map = Gadfly.with_theme(:dark) do
    Gadfly.plot(cells_meanDF[(~).(isnan.(cells_meanDF.df_f_mean)), :], x=:x, y=:y, Gadfly.Geom.point, size=[targetSizePx], color=:df_f_mean,
        Gadfly.Scale.color_continuous(minvalue=-1, maxvalue=1),
        Gadfly.Coord.cartesian(yflip=true, fixed=true))
end

# img = PNG(joinpath(plotDir, "$(expName)_df_f_map.png"), 6inch, 5inch)
img = SVG(joinpath(plotDir, "$(expName)_df_f_map.svg"), 6inch, 5inch)
Gadfly.draw(img, p_df_map)
p_df_map

## extract and save traces for easy plotting later...
fluor = DataFrame(time=Float64[], f=Float64[], cellID=UInt32[], stimStart=UInt32[],
    laserPower=Float64[], stimFreq=Float64[])

# TODO also compare on target vs off-target traces...? Number of off-target change
# as function of # of cell stim?

nTime = before + maximum(cells.stimStop - cells.stimStart) + after
# TODO add @threads ?
for cell in eachrow(cells)
    x, y, z, stimStart, stimStop, laserPower, stimFreq, cellID = cell[[:x, :y, :z, :stimStart, :stimStop, :laserPower, :stimFreq, :cellID]]
    roiM = roiMask[cellID]
    theEnd = minimum([stimStop+after, size(tseries, ndims(tseries))])
    theStart = maximum([stimStart-before, 1])
    plotRange = theStart:theEnd
    timeRange = (plotRange .- stimStart) ./ volRate
    fluorescentTrace = extractTrace(tseries[:,:,:,plotRange], roiM)
    for (t,f) in zip(timeRange, fluorescentTrace)
        push!(fluor, (time=t, f=f, cellID=cellID, stimStart=stimStart,
            stimFreq=stimFreq, laserPower=laserPower ./ 1mW))
    end
end
first(fluor,5)

fluor_arrow_path = tseriesDir*"_cells_fluorescence.arrow"
if ~isfile(fluor_arrow_path)
    open(fluor_arrow_path, "w") do io
        Arrow.write(io, fluor)
    end
end

## save red
red_tseries = loadTseries(tseriesDir,"Ch2");
##
avg_red = mean(red_tseries[:,:,:,1:500],dims=4)[:,:,:,1];
red = RGB.(imadjustintensity(adjust_gamma(avg_red[:,:,5], 0.5)))
channelview(red)[[2,3],:,:] .= 0
red


## pick cells to stim
maxT = minimum([5000, size(tseries,4)])
avg_vol = mean(tseries[:,:,:,1:10:maxT], dims=4)[:,:,:,1];

z = 10
stimNum = 2
f = mean(avgStim[:,:,:,stimNum,end-window+1:end],dims=4)[:,:,:,1]
f0 = mean(avgStim[:,:,:,stimNum,1:window],dims=4)[:,:,:,1]
df = f - f0
df_f = df./f0;

using AbstractPlotting, GLMakie
## we manually rerun this cell for each plane & save the clicks into
# an array in the next cell
# img = df_f[:,:,z]
img = avg_vol[:,:,10]
img .= img[end:-1:1,:]
img = permutedims(img, [2,1])
scene = Scene(resolution = size(img))
AbstractPlotting.heatmap!(scene, img,
    # scale_plot = false,
)
    # colorrange=(-2.5,2.5),colormap=:vik)
clicks = Node(Point2f0[(0,0)])
on(scene.events.mousebuttons) do buttons
   if ispressed(scene, Mouse.left)
       pos = to_world(scene, Point2f0(scene.events.mouseposition[]))
       clicks[] = push!(clicks[], pos)
   end
   return
end
AbstractPlotting.scatter!(scene, clicks, color = :green, marker = '+', markersize = 10)

# Do not execute beyond this point!

RecordEvents(scene, "output")

##
z4_points_s1 = clicks[][2:end]
z4_points_s2 = clicks[][2:end]
z5_points_s1 = clicks[][2:end]
z6_points_s1 = clicks[][2:end-1]
z8_points_s2 = clicks[][2:end]
z9_points_s2 = clicks[][2:end]
z10_points_s2 = clicks[][2:end]
raphe_points = clicks[][2:end]
##

ipn_targets = aa2a(collect.(vcat(z4_points_s1, z4_points_s2, z5_points_s1, z6_points_s1,
    z8_points_s2, z9_points_s2)))
ipn_z = zeros(size(ipn_targets,1))

start = 1
n = (size(z4_points_s1,1)+size(z4_points_s2,1))
ipn_z[start:start+n-1] .= etlVals[4] * 10^-6
start += n
n = size(z5_points_s1,1)
ipn_z[start:start+n-1] .= etlVals[5] * 10^-6
start += n
n = size(z6_points_s1,1)
ipn_z[start:start+n-1] .= etlVals[6] * 10^-6
start += n
n = size(z8_points_s2,1)
ipn_z[start:start+n-1] .= etlVals[8] * 10^-6
start += n
n = size(z9_points_s2,1)
ipn_z[start:start+n-1] .= etlVals[9] * 10^-6
start += n
@assert start - 1 == size(ipn_targets,1)

# account for Makie weirdness
ipn_targets = ipn_targets[:,[2,1]]
ipn_targets[:,1] = 512 .- ipn_targets[:,1]

ipn_targets = hcat(ipn_targets,ipn_z)

raphe_targets = aa2a(collect.(raphe_points))
# account for Makie weirdness
raphe_targets = raphe_targets[:,[2,1]]
raphe_targets[:,1] = 512 .- raphe_targets[:,1]
raphe_z = etlVals[10] * 10^-6
raphe_targets = hcat(raphe_targets,fill(raphe_z, size(raphe_targets,1)))
##

target_groups = [ipn_targets[:,[2,1,3]], raphe_targets[:,[2,1,3]]]
@show size.(target_groups)

N = sum(size.(target_groups,1))
num_groups = length(target_groups)
## Save files for SLM stim
powers = [1]
nPowers = length(powers)
frequencies = repeat([5], nPowers)
name_str = "IPN_Raphe"
name = "$(num_groups)groups_$(N)cells_$name_str"
expSLMdir = joinpath(fishDir,"slm")
outname = joinpath(expSLMdir, name)
if ~isdir(expSLMdir)
    mkdir(expSLMdir)
end
if isfile(outname*".txt")
    # error("file already exists! refusing to clobber")
end

# create_slm_stim(target_groups, outname,
#     localToRemote = matpath -> "Y:" * replace(matpath[15:end], "/" => "\\"),
#     powers=collect(1:nPower)/nPower)

create_slm_stim(target_groups, outname,
    # 9 is path "/scratch/" includes trailing
    localToRemote = matpath -> "T:" * replace(matpath[9:end], "/" => "\\"),
    powers=powers, frequencies=frequencies, slmNum=slmNum)

##
locs = aa2a(collect.(z5_points_s1))
locs = locs[:,[2,1]]
locs[:,1] = 512 .- locs[:,1]

tempavgImageAdj = adjust_gamma(imadjustintensity(avg_vol[:,:,10]), 0.5)
tempavgImageAdj = RGB.(tempavgImageAdj)
channelview(tempavgImageAdj)[[1,3],:,:,:] .= 0

avgImgWithTargets = addTargetsToImage(copy(tempavgImageAdj),
    # Int.(round.(locs,digits=0)),
    Int.(round.(vcat(target_groups...)[:,[2,1]],digits=0)),
    # Int.(round.(raphe_targets,digits=0)),
    targetSize=targetSizePx)