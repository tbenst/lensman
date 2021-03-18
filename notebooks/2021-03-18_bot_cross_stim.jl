## analyze cross-stim from imaging laser
# ENV["DISPLAY"] = "localhost:12"
using Lensman, ImageView, Glob, LibExpat, StatsBase, ProgressMeter, DataFrames,
    MAT, Images, ImageDraw, Arrow, CSV
import Base.Threads.@threads
using Unitful: μm, m, s
import Plots
using Gadfly
import PyPlot
plt = PyPlot
matplotlib = plt.matplotlib

slmDir = "/mnt/b115_mSLM/mSLM/SetupFiles/Experiment/"
# matched...
gcampDir = "/scratch/b115/2021-03-16_h2b6s/fish1"
h33rDir = "/scratch/b115/2021-03-16_h33r-chrmine_h2b6s/fish4"
rsDir = "/scratch/b115/2021-03-16_rschrmine_h2b6s/fish3/"
wtDir = "/scratch/b115/2021-03-16_wt-chrmine_h2b6s/fish2/"

genotypes = ["Control", "H33R", "rsChRmine", "WT"]
dirs = [gcampDir, h33rDir, rsDir, wtDir]

for (g,d) in zip(genotypes, dirs)
    tifDir = glob("*100x100*", d)[end]
    tifDir = glob("*512x512*", d)[end]
    # tifDir = glob("*10plane*", d)[end] # missing BOT for H33R
    csvs = glob("*botData.csv", tifDir)
    global bot
    global df
    df = nothing
    for (t,csv) in enumerate(csvs)
        bot = CSV.File(open(read, csv)) |> DataFrame
        insertcols!(bot, "trial_num" => t, "genotype" => g)
        if isnothing(df)
            df = bot
        else
            df = vcat(df, bot)
        end
    end
    df
end
df
##
savePrefix = "100x100_forebrain_cross-stim_"

plotDir = "/home/tyler/Dropbox/dlab/plots/2021-03-16"
if ~isdir(plotDir)
    mkdir(plotDir)
end

##
readdir(gcampDir)

##
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
    # TODO: mean, ymin, ymax are duplicated...?
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
    for df in groupby(fluor, [:genotype, :trial, :imagingPower])
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

df_f = DataFrame(time=Float64[], imagingPower=String[], mean=Float64[],
    ymin=Float64[], ymax=Float64[], genotype=String[])
    
for df in groupby(rawToDf_f(df), [:genotype, :time, :imagingPower])
    time = df[1,:time]
    imagingPower = df[1,:imagingPower]
    if imagingPower==140
        imagingPower = "10mW"
    elseif imagingPower==125
        imagingPower = "8mW"
    end
    genotype = df[1,:genotype]
    df_fs = df.df_f
    the_mean = mean(df_fs)
    ymin = percentile(df_fs,2.5)
    ymax = percentile(df_fs,97.5)
    push!(df_f, (mean=the_mean, imagingPower=imagingPower, time=time,
        ymin=ymin, ymax=ymax, genotype=genotype))
end
##

line = plot(df_f, x=:time,
    xgroup=:imagingPower,
    y=:mean,
    ymin=:ymin,
    ymax=:ymax,
    color=:genotype,
    Geom.subplot_grid(Geom.line, Geom.ribbon))

Gadfly.draw(PNG(joinpath("/mnt/deissero/users/tyler/plots/","2021-02-23_forebrain_cross-stim-rs-wt-control.png"), 9inch, 5inch), line)
