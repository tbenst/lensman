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

crossStimDF = nothing
genotypes = ["Control", "H33R", "rsChRmine", "WT"]
dirs = [gcampDir, h33rDir, rsDir, wtDir]

for (g,d) in zip(genotypes, dirs)
    tifDir100 = glob("*100x100*", d)[end]
    tifDir512 = glob("*512x512*", d)[end]
    @show tifDir100, tifDir512
    # tifDir = glob("*10plane*", d)[end] # missing BOT for H33R
    csvs100 = glob("*botData.csv", tifDir100)
    csvs512 = glob("*botData.csv", tifDir512)
    global bot
    global df
    for (path,r) in zip([csvs100, csvs512], [100, 512])
        for (t,csv) in enumerate(path)
            bot = CSV.File(open(read, csv)) |> DataFrame
            if r == 512
                # bot = bot[7:10:end,:] # only plane 7 is our ROI
                # rsChRmine 512: ROI2
                # H2B control 512: ROI2
                # H33R 512: ROI2
                # WT 512: ROI2
                bot[:,"Region 1"] .= bot[:,"Region 2"] # hack to use same left habenula ROI as roi=100
            end
            insertcols!(bot, "trial_num" => t, "genotype" => g, "roi" => r, "time" => 1:size(bot,1)*1.0)

            if isnothing(crossStimDF)
                crossStimDF = bot
            else
                if length(names(bot)) != length(names(crossStimDF))
                    for n in unique([names(bot); names(crossStimDF)]), df in [bot, crossStimDF]
                        n in names(df) || (df[!,n] .= missing)
                    end
                end
                crossStimDF = vcat(crossStimDF, bot)
            end
        end
    end
end


crossStimDF

@warn "hack hardcode to fix time of 20s"
idx100 = crossStimDF.roi .== 100
idx512 = crossStimDF.roi .== 512
crossStimDF.time[idx100] ./= (maximum(crossStimDF.time[idx100])/20)
crossStimDF.time[idx512] ./= (maximum(crossStimDF.time[idx512])/20)

##

plotDir = "/home/tyler/Dropbox/dlab/plots/2021-03-16"
if ~isdir(plotDir)
    mkdir(plotDir)
end
##

meanDF = combine(groupby(crossStimDF[:, ["time", "Region 1", "trial_num", "genotype", "roi"]],
    ["time", "genotype", "roi"]), "Region 1" => mean, "Region 1" => maximum, "Region 1" => minimum)

##
line = plot(meanDF[meanDF.roi .== 100, :], x=:time,
    xgroup=:genotype,
    y="Region 1_mean",
    ymin="Region 1_minimum",
    ymax="Region 1_maximum",
    Geom.subplot_grid(Geom.line, Geom.ribbon))
    # Geom.subplot_grid(Geom.line))


    
Gadfly.draw(SVG(joinpath(plotDir,"cross-stim-raw.svg"), 9inch, 5inch), line)
# Gadfly.draw(PNG(joinpath(plotDir,"$savePrefix.png"), 9inch, 5inch), line)
line

##

function rawToDf_f(fluor::DataFrame; before=50)
    fluorDF = []
    # operate over each cell's trace...
    for df in groupby(fluor, [:genotype, :trial_num, :roi])
        kalmanFilt = imageJkalmanFilter(df[:, "Region 1"])
        f0 = mean(kalmanFilt[1:before])
        df = copy(df)
        df_f = @. (kalmanFilt - f0) / (f0 + 10)
        insertcols!(df, size(df,2), :df_f => df_f)
        select!(df, Not("Region 1"))
        push!(fluorDF, df)
    end
    fluorDF = vcat(fluorDF...)
end

cross_stim_df_f = rawToDf_f(crossStimDF, before=50)

mean_cross_stim_df_f = combine(groupby(cross_stim_df_f[:, ["time", "df_f", "trial_num", "genotype", "roi"]],
    ["time", "genotype", "roi"]), "df_f" => mean, "df_f" => maximum, "df_f" => minimum)


line = plot(mean_cross_stim_df_f[meanDF.roi .== 100, :], x=:time,
    xgroup=:genotype,
    y="df_f_mean",
    ymin="df_f_minimum",
    ymax="df_f_maximum",
    Geom.subplot_grid(Geom.line, Geom.ribbon))
    # Geom.subplot_grid(Geom.line))


    
Gadfly.draw(SVG(joinpath(plotDir,"cross-stim-df_f-100.svg"), 9inch, 5inch), line)
# Gadfly.draw(PNG(joinpath(plotDir,"$savePrefix.png"), 9inch, 5inch), line)
line

##
cross_stim_df_f = rawToDf_f(crossStimDF, before=5)

mean_cross_stim_df_f = combine(groupby(cross_stim_df_f[:, ["time", "df_f", "trial_num", "genotype", "roi"]],
    ["time", "genotype", "roi"]), "df_f" => mean, "df_f" => maximum, "df_f" => minimum)


line = plot(mean_cross_stim_df_f[meanDF.roi .== 512, :], x=:time,
    xgroup=:genotype,
    y="df_f_mean",
    ymin="df_f_minimum",
    ymax="df_f_maximum",
    Geom.subplot_grid(Geom.line, Geom.ribbon))
    # Geom.subplot_grid(Geom.line))


    
Gadfly.draw(SVG(joinpath(plotDir,"cross-stim-df_f-512.svg"), 9inch, 5inch), line)
# Gadfly.draw(PNG(joinpath(plotDir,"$savePrefix.png"), 9inch, 5inch), line)
line