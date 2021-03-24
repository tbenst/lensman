## analyze cross-stim from imaging laser
# ENV["DISPLAY"] = "localhost:12"
using Lensman, ImageView, Glob, LibExpat, StatsBase, ProgressMeter, DataFrames,
    MAT, Images, ImageDraw, Arrow, CSV, ImageMagick, Images, HypothesisTests,
    DSP
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
genotypes = ["Control", "H33R", "rsChRmine", "ChRmine"]
dirs = [gcampDir, h33rDir, rsDir, wtDir]

@warn "we hardcode plane 7 for 10 plane"

for (g,d) in zip(genotypes, dirs)
    tifDir100 = glob("*100x100*", d)[end]
    tifDir512 = glob("*512x512*", d)[end]
    tifDir5120 = glob("*10plane*", d)[end]
    @show tifDir100, tifDir512, tifDir5120
    # tifDir = glob("*10plane*", d)[end] # missing BOT for H33R
    csvs100 = glob("*botData.csv", tifDir100)
    csvs512 = glob("*botData.csv", tifDir512)
    csvs5120 = glob("*botData.csv", tifDir5120)
    global bot
    global df
    for (path,r) in zip([csvs100, csvs512, csvs5120], [100, 512, 5120])
    # for (path,r) in zip([csvs5120], [5120])
        cur_time_idx = 1 # used for z-series only
        prev_time = -Inf
        cur_trial_idx = 0
        for (t,csv) in enumerate(path)
            bot = CSV.File(open(read, csv)) |> DataFrame
            if (r == 512) | (r == 5120)
                # bot = bot[7:10:end,:] # only plane 7 is our ROI
                # rsChRmine 512: ROI2
                # H2B control 512: ROI2
                # H33R 512: ROI2
                # WT 512: ROI2
                bot[:,"Region 1"] .= bot[:,"Region 2"] # hack to use same left habenula ROI as roi=100
            end
            if r == 5120 # volumetric imaging of 10 planes
                # for zseries, each .csv is 1 z-stack
                if bot.Timestamp[1] - prev_time > 1
                    cur_time_idx = 1
                    cur_trial_idx += 1
                end
                prev_time = bot.Timestamp[end] # for next loop
                bot = bot[[7],:]
                insertcols!(bot, "trial_num" => cur_trial_idx, "genotype" => g, "roi" => r,
                    "time" => cur_time_idx)
                cur_time_idx += 1
            else
                insertcols!(bot, "trial_num" => t, "genotype" => g, "roi" => r,
                    "time" => 1:size(bot,1)*1.0)
            end

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
##
@warn "hack hardcode to fix time to 20s"
idx100 = crossStimDF.roi .== 100
idx512 = crossStimDF.roi .== 512
idx5120 = crossStimDF.roi .== 5120
crossStimDF.time[idx100] ./= (maximum(crossStimDF.time[idx100])/20)
crossStimDF.time[idx512] ./= (maximum(crossStimDF.time[idx512])/20)
crossStimDF.time[idx5120] ./= (maximum(crossStimDF.time[idx5120])/20)

##

plotDir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure"
if ~isdir(plotDir)
    mkdir(plotDir)
end
##

function max_95ci(x)
    confint(OneSampleTTest(x))[2]
end

function min_95ci(x)
    confint(OneSampleTTest(x))[1]
end

meanDF = combine(groupby(crossStimDF[:, ["time", "Region 1", "trial_num", "genotype", "roi"]],
    ["time", "genotype", "roi"]), "Region 1" => mean, "Region 1" => max_95ci, "Region 1" => min_95ci)

##

chrmine_paper = Theme(
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


line = plot(meanDF[meanDF.roi .== 100, :], x=:time,
    xgroup=:genotype,
    y="Region 1_mean",
    # ymin="Region 1_minimum",
    # ymax="Region 1_maximum",
    ymin="Region 1_min_95ci",
    ymax="Region 1_max_95ci",
    chrmine_paper,
    Geom.subplot_grid(Geom.line, Geom.ribbon))
    # Geom.subplot_grid(Geom.line))


    
Gadfly.draw(SVG(joinpath(plotDir,"cross-stim-raw.svg"), 6cm, 2cm), line)
# Gadfly.draw(PNG(joinpath(plotDir,"$savePrefix.png"), 9inch, 5inch), line)
line

##
"Avoid ringing."
function resample_pad(x, rate)
    # we will drop 10 samples after resample
    reps = Int(ceil(1/rate))*10
    new = vcat(repeat([x[1]],reps), x, repeat([x[end]],reps))
    resample(new,rate)[11:end-10]
end

function rawToDf_f(fluor::DataFrame; before=Dict(100 => 50, 512 => 50, 5120 => 3))
    fluorDF = []
    # operate over each cell's trace...
    for df in groupby(fluor, [:genotype, :trial_num, :roi])
        kalmanFilt = imageJkalmanFilter(df[:, "Region 1"])
        f0 = mean(kalmanFilt[1:before[df.roi[1]]])
        df = copy(df)
        # avoid too much resolution on SVG: resample to 6Hz
        if df.roi[1] == 100
            df = copy(df[1:25:end,:])
            kalmanFilt = resample_pad(kalmanFilt, 0.04)
        elseif df.roi[1] == 512
            df = copy(df[1:5:end,:])
            kalmanFilt = resample_pad(kalmanFilt, 0.2)
        end
        df_f = @. (kalmanFilt - f0) / (f0 + 10)
        insertcols!(df, size(df,2), :df_f => df_f)
        select!(df, Not("Region 1"))
        push!(fluorDF, df)
    end
    fluorDF = vcat(fluorDF...)
end

cross_stim_df_f = rawToDf_f(crossStimDF)

mean_cross_stim_df_f = combine(groupby(cross_stim_df_f[:, ["time", "df_f", "trial_num", "genotype", "roi"]],
    ["time", "genotype", "roi"]), "df_f" => mean, "df_f" => max_95ci, "df_f" => min_95ci)

##
sz = 100
sz = 512
# sz = 5120
df = mean_cross_stim_df_f[mean_cross_stim_df_f.roi .== sz, :]
df = df[df.genotype .!= "H33R", :]
line = plot(df, x=:time,
    xgroup=:genotype,
    y="df_f_mean",
    # ymin="df_f_minimum",
    # ymax="df_f_maximum",
    ymin="df_f_min_95ci",
    ymax="df_f_max_95ci",
    Guide.ylabel("df/f"),
    chrmine_paper,
    Geom.subplot_grid(Geom.line, Geom.ribbon, Gadfly.Coord.Cartesian(ymin=-0.1, ymax=0.5)))
    # Geom.subplot_grid(Geom.line))


plot_path = joinpath(plotDir,"cross-stim-df_f-$sz.svg")
Gadfly.draw(SVG(plot_path, 13.5cm, 6cm), line)
@show plot_path
# Gadfly.draw(PNG(joinpath(plotDir,"$savePrefix.png"), 12cm, 4cm), line)
line

##
# compare expression
function get_chrmine_intensity(fishDir, k=x->Float64(mean(x[200:300,200:300])))
    tif = glob("*Ch2*", glob("*115*new", fishDir)[end])[1]
    tif = reinterpret(UInt16,ImageMagick.load(tif))
    k(tif)
end

# collect(zip(dirs, map(get_chrmine_intensity, dirs)))
# collect(zip(dirs, map(get_chrmine_intensity, dirs)))
# collect(zip(dirs, map(x -> get_chrmine_intensity(x,
#     x->Float64(mean(x[160:310, 160:310]))), dirs))) 
collect(zip(dirs, map(x -> get_chrmine_intensity(x,
    x->Float64(mean(x[370:380, 278:288]))), dirs))) 
#  ("/scratch/b115/2021-03-16_h2b6s/fish1", 2.5041322314049586)
#  ("/scratch/b115/2021-03-16_h33r-chrmine_h2b6s/fish4", 5.867768595041323)
#  ("/scratch/b115/2021-03-16_rschrmine_h2b6s/fish3/", 371.64462809917353)
#  ("/scratch/b115/2021-03-16_wt-chrmine_h2b6s/fish2/", 157.26446280991735)
## average fluorescence in left habenula
i = 3
tif = glob("*Ch2*", glob("*115*new", dirs[i])[end])[1]
tif = reinterpret(UInt16,ImageMagick.load(tif))
avg = Float64(mean(tif[178:178+25,431:431+25]))
@show genotypes[i], avg
i = 4
tif = glob("*Ch2*", glob("*115*new", dirs[i])[end])[1]
tif = reinterpret(UInt16,ImageMagick.load(tif))
avg = Float64(mean(tif[226:251,442:442+25])) * 1.6
@show genotypes[i], avg
i = 1
tif = glob("*Ch2*", glob("*115*new", dirs[i])[end])[1]
tif = reinterpret(UInt16,ImageMagick.load(tif))
avg = Float64(mean(tif[226:251,442:442+25]))
@show genotypes[i], avg


##
# safe values:
maximum(map(x -> get_chrmine_intensity(x,
    y->Float64(percentile(adjust_gamma(y, 2.)[:],99))), dirs))

    # y->percentile(y[:],95)), dirs)))
##
for (i,g) in enumerate(genotypes)
    d = dirs[i]
    tif = glob("*Ch2*", glob("*115*new", d)[end])[1]
    tif = ImageMagick.load(tif)
    # adjust_histogram(tif, LinearStretching(dst_minval = 0, dst_maxval = 1))
    @warn "hardcode fix for different pockels value (eye brightness should be the same)."
    if i == 4
        # approximately (800/750)^7.5
        channelview(tif) .*= 1.6
    else
    end
    tif = RGB.(clamp01.(adjust_gamma(tif, 2.)*5))
    # @show percentile(channelview(tif)[:],95)
    channelview(tif)[[2,3],:,:] .= 0
    tif_path = joinpath(plotDir,"$g-1115nm.png")
    save(tif_path, tif)
end
# adjust_gamma(imadjustintensity(tif),3)


##
d = dirs[3]
@show d
tif = glob("*Ch2*", glob("*115*new", d)[end])[1]
tif = ImageMagick.load(tif)
# channelview(tif) .*= 2
# mean(channelview(adjust_gamma(tif*one(tif[1])*2, 2.)*3))
# mean(channelview(adjust_gamma(tif, 2.)*3))
tif = RGB.(clamp01.(adjust_gamma(tif, 2.)*3))
tif
# tif = RGB.(clamp01.(adjust_gamma(tif, 2.)*3))
# tif .*= 5
# Float64(RGB.(tif .* 1)