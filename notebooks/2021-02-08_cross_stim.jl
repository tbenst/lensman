## analyze cross-stim from imaging laser
ENV["DISPLAY"] = "localhost:12"
using Lensman, ImageView, Glob, LibExpat, StatsBase, ProgressMeter, DataFrames,
    MAT, Images, ImageDraw, Arrow
import Base.Threads.@threads
using Unitful: μm, m, s
import Plots
using Gadfly
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib

slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"

##
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-02_wt_chrmine_GC6f/fish3/"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-02_f1_h33r_GC6f_6dpf/fish2/"
fishDir = "/mnt/deissero/users/tyler/b115/2021-02-08_gcamp6f_6dpf/fish1"
# "/data/dlab/b115/2021-02-09_wt-chrmine_gcamp6f_7dpf/"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-02_f1_h33r_GC6f_6dpf/fish1_nochrmine"
# tseriesDirs = glob("TSeries_cross-stim*", fishDir)
tseriesDirs = glob("TSeries-cross-stim*", fishDir)


plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end


expName = splitpath(tseriesDirs[1])[end]
expXmlPath = joinpath(tseriesDirs[1],expName*".xml")
expDate, frameRate, etlVals = getExpData(expXmlPath)
targetsDir = joinpath(fishDir, "slm")
matFiles = findMatGroups(targetsDir, group_re=Regex(raw".*(\d+)\.mat$"))
target_groups = [mat["experimentS"]["targets"][1]
    for mat in matread.(matFiles)];

Z = 1
H = 512
W = 512
##
microscope_units_xy = 0.6299544139175637μm
target_size_px = 7μm * (14.4 / 25) / microscope_units_xy
stim_mask = Gray.(zeros(Bool, H, W))
for (stim_idx, target_group) in enumerate(target_groups)
    for (x, y, z) in eachrow(target_group)
        x = round(Int,x)
        y = round(Int,y)
        draw!(view(stim_mask, :, :), Ellipse(CirclePointRadius(x, y, target_size_px)))
    end
end
stim_mask = convert(Array{Bool}, stim_mask)
@show sum(stim_mask)
Gray.(stim_mask)



##
df = DataFrame(time=Float64[], f=Float64[], imagingPower=Int64[], trial=Int64[],
               mean=Float64[], ymin=Float64[], ymax=Float64[])
@showprogress for tseriesDir in tseriesDirs
    # tseriesDir = tseriesDirs[1]
    expName = splitpath(tseriesDir)[end]
    expXmlPath = joinpath(tseriesDir,expName*".xml")

    tseries_xml = parseXML(expXmlPath);
    expDate, frameRate, etlVals = getExpData(expXmlPath)
    imagingPockels = round(Int,getImagingPockels(tseries_xml))

    tseries = loadTseries(tseriesDir);
    H, W, T, nTrial = size(tseries)

    # use only the stimulation mask
    # avgTrace = reshape(extractTrace(reshape(tseries,H,W,T*nTrial), stim_mask), T, nTrial)
    avgTrace = mean(tseries, dims=(1,2))[1,1,:,:]
    ts = collect(1:T) / frameRate
    for (t,trace) in enumerate(eachcol(avgTrace))
        df = vcat(df,
            DataFrame(time=ts, f=trace, imagingPower=imagingPockels,
                trial=repeat([t],length(trace)), mean=mean(avgTrace,dims=2)[:,1],
                    ymin=percentile.(eachrow(avgTrace),2.5),
                    ymax=percentile.(eachrow(avgTrace),97.5)))
    end
end
##
line = plot(df, x=:time,
    xgroup=:imagingPower,
    y=:mean,
    ymin=:ymin,
    ymax=:ymax,
    Geom.subplot_grid(Geom.line, Geom.ribbon))

Gadfly.draw(SVG(joinpath(plotDir,"cross_stim.svg"), 9inch, 5inch), line)
Gadfly.draw(PNG(joinpath(plotDir,"cross_stim.png"), 9inch, 5inch), line)
line
##
tseriesDir = tseriesDirs[1]
expName = splitpath(tseriesDir)[end]
tseries = loadTseries(tseriesDir);
expXmlPath = joinpath(tseriesDir,expName*".xml")

tseries_xml = parseXML(expXmlPath);
expDate, frameRate, etlVals = getExpData(expXmlPath)
imagingPockels = round(Int,getImagingPockels(tseries_xml))



# in ChrMine habenula lights up but that's it...?
# but also in the control fish...? (but was injected...)
trialAvg = mean(tseries, dims=(4))[:,:,:,1]
imshow(trialAvg)

##
f0 = mean(trialAvg[:,:,1:10],dims=3)[:,:,1] .+ 1
df_fs = []
for i = 1:3
    s = i*10+1
    e = (i+1)*10
    f = mean(trialAvg[:,:,s:e],dims=3)[:,:,1]
    df = f - f0
    df_f = df./f0
    push!(df_fs, df_f)
end

f = mean(trialAvg[:,:,end-10:end],dims=3)[:,:,1]
df = f - f0
df_f = df./f0



cmax = 2
cmin = -0.75
cnorm = matplotlib.colors.TwoSlopeNorm(vmin=cmin,vcenter=0,vmax=cmax)
# fig = plt.figure(figsize=(30,10))
fig = plt.figure()
plt.axis("off")
# plt.imshow(hcat(df_fs...), cmap="RdBu_r", norm=cnorm)
plt.imshow(df_f, cmap="RdBu_r", norm=cnorm)
plt.colorbar()
fig.savefig(joinpath(plotDir,"max_power_cross_stim_df_map_p$imagingPockels.png"))
joinpath(plotDir,"max_power_cross_stim_df_map_p$imagingPockels.png")

##
open(joinpath(fishDir, "cross-stim_df.arrow"), "w") do io
    Arrow.write(io, df)
end
