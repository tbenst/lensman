## analyze cross-stim from imaging laser
# ENV["DISPLAY"] = "localhost:12"
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
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-22_h2b6s/fish2/" # habenula shows cross stimulation
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-22_h2b6s_wt-chrmine/fish1/" # most are hab, but 4th may be forebrain...?

fishDir = "/mnt/deissero/users/tyler/b115/2021-02-22_h2b6s_wt-chrmine/fish3/"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2"
# fishDir = "/mnt/deissero/users/tyler/b115/2021-02-23_h2b6s/fish1"

tseriesDirs = glob("TSeries-cross-stim*", fishDir)
# tseriesDirs = glob("TSeries-cross-stim*", fishDir)[end-1:end]

savePrefix = "100x100_forebrain_cross-stim_"
# tseriesDirs = tseriesDirs[[end]]

plotDir = joinpath(fishDir, "plots")
if ~isdir(plotDir)
    mkdir(plotDir)
end


expName = splitpath(tseriesDirs[1])[end]
expXmlPath = joinpath(tseriesDirs[1],expName*".xml")
expDate, frameRate, etlVals = getExpData(expXmlPath)

##
df = DataFrame(time=Float64[], f=Float64[], imagingPower=Int64[], trial=Int64[],
               mean=Float64[], ymin=Float64[], ymax=Float64[])
avgImages = []
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
    push!(avgImages, mean(tseries, dims=(3,4))[:,:,1,1])
    ts = collect(1:T) / frameRate
    # TODO: mean, ymin, ymax seem wrong / duplicate...?
    for (t,trace) in enumerate(eachcol(avgTrace))
        df = vcat(df,
            DataFrame(time=ts, f=trace, imagingPower=imagingPockels,
                trial=repeat([t],length(trace)), mean=mean(avgTrace,dims=2)[:,1],
                    ymin=percentile.(eachrow(avgTrace),2.5),
                    ymax=percentile.(eachrow(avgTrace),97.5)))
    end
end

open(joinpath(fishDir, "$(savePrefix)_fluorescence.arrow"), "w") do io
    Arrow.write(io, df)
end

# TODO: these images seem subtly wrong...?!? some weird mirroring ?!?!?
im = RGB.(imadjustintensity(adjust_gamma(avgImages[1], 0.5)))
channelview(im)[[1,3],:,:] .= 0
save(joinpath(fishDir, savePrefix*".png"), im)
im


##
line = plot(df, x=:time,
    xgroup=:imagingPower,
    y=:mean,
    ymin=:ymin,
    ymax=:ymax,
    Geom.subplot_grid(Geom.line, Geom.ribbon))

Gadfly.draw(SVG(joinpath(plotDir,"$savePrefix.svg"), 9inch, 5inch), line)
Gadfly.draw(PNG(joinpath(plotDir,"$savePrefix.png"), 9inch, 5inch), line)
line

##

arrowWT = glob("*.arrow", "/mnt/deissero/users/tyler/b115/2021-02-22_h2b6s_wt-chrmine/fish3/")[1]
arrowRS = glob("*.arrow", "/mnt/deissero/users/tyler/b115/2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2")[1]
arrowControl = glob("*.arrow", "/mnt/deissero/users/tyler/b115/2021-02-23_h2b6s/fish1")[1]

dfWT = copy(Arrow.Table(arrowWT) |> DataFrame)
insertcols!(dfWT, size(dfWT,2)+1, :genotype => "WT")
dfRS = copy(Arrow.Table(arrowRS) |> DataFrame)
insertcols!(dfRS, size(dfRS,2)+1, :genotype => "RS")
dfControl = copy(Arrow.Table(arrowControl) |> DataFrame)
insertcols!(dfControl, size(dfControl,2)+1, :genotype => "Control")
df = vcat(dfWT, dfRS, dfControl)
##
function rawToDf_f(fluor::DataFrame, before=50)
    fluorDF = []
    # operate over each cell's trace...
    for df in groupby(fluor, [:genotype, :trial])
        kalmanFilt = imageJkalmanFilter(df.f)
        f0 = mean(kalmanFilt[1:before])
        df = copy(df)
        df_f = @. (kalmanFilt - f0) / (f0 + 10)
        insertcols!(df, size(df,2), :df_f => df_f)
        select!(df, Not(:f))
        push!(fluorDF, df)
    end
    fluorDF = vcat(fluorDF...)
end

df_f = rawToDf_f(df)
insertcols!(dfControl, size(dfControl,2)+1, :genotype => "Control")
##

line = plot(df_f, x=:time,
    xgroup=:imagingPower,
    y=:df_f,
    # ymin=:ymin,
    # ymax=:ymax,
    Geom.subplot_grid(Geom.line, Geom.ribbon))
