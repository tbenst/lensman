# ENV["DISPLAY"] = "localhost:11"
# mostly same as 2021-03-04 notebook, but updated for latest
# DataFrames
##
using AlgebraOfGraphics, CairoMakie
import AbstractPlotting
using Makie.GeometryBasics
# import Cairo, Fontconfig
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Dates, ImageDraw, MAT, StatsBase,
    ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions, Fontconfig
# import Gadfly
using Unitful: μm, m, s
aog = AlgebraOfGraphics
# PNG

## compare ChRmine: WT, RS & H33R

# tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
tseriesRootDir = "/mnt/deissero/users/tyler/b115"
# tseriesRootDir = "/scratch/b115"
# h33r1 = "$tseriesRootDir/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-256cell-8concurrent-4freq-055"
# h33r1 = h33r1*"_cells_fluorescence.arrow"
# wt_path = "$tseriesRootDir/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043"
# tseriesDir = joinpath(tseriesRootDir, "2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2/TSeries-128cell-4concurrent-3power-skip7-044")


# round 1
# h33r_path = joinpath(tseriesRootDir, "2021-02-16_6f_h33r_f0_6dpf/fish2/TSeries-256cell-8concurrent-4freq-051")
# rs_path = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043")
# rs_path = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish2/TSeries-1024cell-32concurrent-048")
# wt_path = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-10zplane-077"
# wt_path = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-076"

# round 2
h33r_path = "$tseriesRootDir/2021-03-16_h33r-chrmine_h2b6s/fish4/TSeries_64cell_8concurrent_2power_8rep-607"
rs_path = "$tseriesRootDir/2021-03-16_rschrmine_h2b6s/fish3/TSeries_64cell_8concurrent_2power_8rep-407"
gcamp_path = "$tseriesRootDir/2021-03-16_h2b6s/fish1/TSeries_64cell_8concurrent_2power_8rep-207"
wt_path = "$tseriesRootDir/2021-03-16_wt-chrmine_h2b6s/fish2/TSeries_64cell_8concurrent_2power_8rep-221"
# check pulse rate
resources = Resources();
rwt = Recordings[
    "2021-03-16_wt-chrmine_h2b6s/fish2/TSeries_64cell_8concurrent_2power_8rep-221"
](;resources...,tseries_read_strategy=:lazy_tiff);
@show rwt[:nstim_pulses] # 10
#

function readExp(tseries_dir, genotype, cellIDnonce)
    cells = tseries_dir*"_cellsDF.arrow"
    fluor = tseries_dir*"_cells_fluorescence.arrow"
    cells = copy(Arrow.Table(cells) |> DataFrame);
    fluor = copy(Arrow.Table(fluor) |> DataFrame);
    insertcols!(fluor, 1+size(fluor,2), :genotype => genotype)
    if ~("stimFreq" in names(fluor))
        insertcols!(fluor, :stimFreq =>
        [cells[findall(cells.stimStart .== start)[1],:stimFreq]
        for start in fluor.stimStart])
    end
    if ~("laserPower" in names(fluor))
        insertcols!(fluor, :laserPower =>
        [cells[findall(cells.stimStart .== start)[1],:laserPower]
            for start in fluor.stimStart])
    end
    cells.cellID .+= cellIDnonce
    fluor.cellID .+= cellIDnonce
    cells, fluor
end

wt_cells, wt = readExp(wt_path, "ChRmine", 0e4)
h33r_cells, h33r = readExp(h33r_path, "fChRmine", 1e4)
rs_cells, rs = readExp(rs_path, "rsChRmine", 2e4)
# gcamp_cells, gcamp = readExp(gcamp_path, "control", 3e4)


# fluor = vcat(h33r, rs, wt, gcamp);
fluor = vcat(h33r, rs, wt);
# plotDir = "$tseriesRootDir/../plots/2021-03-16"
# plotDir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/v3"
plot_dir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/2021-08-27_panels"

##

# downsample data for same sampling rate & smoothing
nTimes = combine(x-> length(unique(x.time)), groupby(fluor, :genotype)).x1
maxTime = combine(x-> maximum(x.time), groupby(fluor, :genotype)).x1
minTime = combine(x-> minimum(x.time), groupby(fluor, :genotype)).x1
frameRates = round.(nTimes ./ (maxTime - minTime), digits=0)
@show frameRates
# resample = Int.(frameRates ./ 3)

# round 1
# resample = [1, 2, 1]
# @warn "resampling hardcode stuff"
# round 2
# resample = [1, 1, 1, 1]
# resample = [1, 2,2,2]
resample = [3, 6,6,6]
# resample = [15, 30,30,30]

genotype_order = map(x->x[1],keys(groupby(fluor, :genotype)))


# length is either 472, 944, or 96....
function rawToDf_f(fluor::DataFrame, genotype_order, resample)
    fluorDF = []
    # operate over each cell's trace...
    @warn "hack hardcoded to adjust timescale"
    for df in groupby(fluor, [:cellID, :stimStart])
        g = df.genotype[1]
        r = resample[findall(genotype_order.==g)[1]]
        if r > 1
            f = df.f[1:r:end-1] .+ df.f[2:r:end]
            # f = df.f[1:2:942] .+ df.f[2:2:942]
            # df = copy(df[1:r:942,:])
            df = copy(df[1:r:end-1,:])
        else
            f = df.f
            df = copy(df[1:r:end,:])
        end
        # @show length(df.f), g, r
        before = 30 # 2 sec
        # kalmanFilt = imageJkalmanFilter(medianfilt(f,3))
        kalmanFilt = imageJkalmanFilter(f)
        f0 = mean(kalmanFilt[1:before])
        df_f = @. (kalmanFilt - f0) / (f0 + 10)
        insertcols!(df, size(df,2), :df_f => df_f)
        select!(df, Not(:f))
        push!(fluorDF, df)
    end
    fluorDF = vcat(fluorDF...)
end

# function rawToDf_f(fluor::DataFrame; before=50)
#     fluorDF = []
#     # operate over each cell's trace...
#     for df in groupby(fluor, [:genotype, :trial_num, :roi])
#         kalmanFilt = imageJkalmanFilter(df[:, "Region 1"])
#         f0 = mean(kalmanFilt[1:before])
#         df = copy(df)
#         df_f = @. (kalmanFilt - f0) / (f0 + 10)
#         insertcols!(df, size(df,2), :df_f => df_f)
#         select!(df, Not(:f))
#         push!(fluorDF, df)
#     end
#     fluorDF = vcat(fluorDF...)
# end


fluorDF = rawToDf_f(fluor, genotype_order, resample)
meanTimeDf = combine(groupby(fluorDF, [:cellID, :time, :stimFreq, :laserPower, :genotype]),
    :df_f => mean, :df_f => minimum, :df_f => maximum)

##
df_fs = DataFrame(cellID=UInt32[], stimFreq=Int64[], laserPower=Float64[],
    genotype=String[], df_f=Float64[])
for df in groupby(meanTimeDf, [:laserPower, :stimFreq, :cellID, :genotype])
    cellID, stimFreq, laserPower, genotype = df[1,[:cellID, :stimFreq, :laserPower, :genotype]]
    idxs = (df.time .> 1.1) .& (df.time .<= 3.1)
    push!(df_fs, (cellID=cellID, stimFreq=stimFreq, laserPower=laserPower,
        genotype=genotype, df_f=mean(df.df_f_mean[idxs])))
end

# Gadfly.plot(df_fs, x=:df_f, ygroup=:genotype,
#     Gadfly.Geom.subplot_grid(Gadfly.Geom.histogram()))

## GCaMP6f

## round 1
h33r_good_idxs = UInt32[10078, 10087, 10184, 10089, 10216, 10077, 10083, 10124, 10200, 10120, 10192, 10226]
# smaller df
# h33r_good_idxs = UInt32[10103, 10083, 10089, 10120, 10078, 10087, 10192, 10186, 10226, 10181, 10127, 10124]
# wt_good_idxs = UInt32[318, 618, 396, 621, 408, 834, 854, 170, 415, 310, 1018, 853]
# wt_good_idxs = UInt32[169, 854, 1024, 307, 834, 1005, 310, 1015, 1016, 318, 853, 439] # 51
# one plane wt
# wt_good_idxs = UInt32[25, 6, 83, 503, 535, 565, 520, 523, 525, 91, 96, 655]
# rs_good_idxs = UInt32[20361, 20341, 20629, 20249, 20641, 20479, 20905, 20960, 20821, 20322, 20373, 20130] #20174
# rs_good_idxs = UInt32[20004, 20087, 20008, 20061, 20025, 20130, 20028, 20149, 20201, 20085, 20135, 20733]


# round 2
# h33r_good_idxs = UInt32[10058, 10005, 10049, 10021, 10047]
# h33r_good_idxs = UInt32[10018, 10044, 10032, 10045, 10046, 10003, 10004, 10060, 10033, 10029, 10012,]
control_good_idxs = UInt32[]
wt_good_idxs = UInt32[42, 14, 9, 47, 49, 12, 32, 43, 13, 61, 41, 57 ]
# 20022, 20055 are weak
rs_good_idxs = UInt32[20009, 20010, 20019, 20004, 20032, 20024, 20064, 20020, 20014, 20017, 20022, 20055]
# 0.5: 20014, 20046, 20017, 20018; 0.3: 20026, 20057; 20020 weird plateau..
# rs_good_idxs = UInt32[20009, 20010, 20024, 20014, 20056, 20064, 20018, 20019, 20046, 20017, 20026, 20057]

# h33r_good_idxs
# wt_good_idxs = wt_good_idxs[1:end-1]
# rs_good_idxs = rs_good_idxs[1:end-1]

@assert length(unique(h33r_good_idxs)) == 12
@assert length(unique(wt_good_idxs)) == 12
@assert length(unique(rs_good_idxs)) == 12
# @assert length(unique(h33r_good_idxs)) == 11
# @assert length(unique(wt_good_idxs)) == 11
# @assert length(unique(rs_good_idxs)) == 11


idxs = vcat(h33r_good_idxs, control_good_idxs, wt_good_idxs, rs_good_idxs)
##
keys(groupby(fluor, [:genotype, :laserPower]))
keys(groupby(meanTimeDf, [:genotype, :laserPower]))
##

# idxs = map(i->i in h33r_good_idxs, meanTimeDf.cellID)
# powerPerGenotype = Dict("ChRmine" => 9.1, "fChRmine" => 16.0,
#     "rsChRmine" => 6.5)
# powerPerGenotype = Dict("ChRmine" => 14.8, "fChRmine" => 16.0,
#     "rsChRmine" => 14.8)
powerPerGenotype = Dict("ChRmine" => 14.8, "fChRmine" => 14.8,
    "rsChRmine" => 14.8)
    # "rsChRmine" => 6.5)
meantraceDF = DataFrame(time=Float64[], mean=Float64[], min=Float64[], max=Float64[], genotype=String[])
for (gen,idxsToUse) in zip(["ChRmine", "rsChRmine"], [wt_good_idxs, rs_good_idxs])
    @show gen
    idxs = map(i->i in idxsToUse, meanTimeDf.cellID)
    idxs = map(i->i in idxsToUse, meanTimeDf.cellID) .&
        (meanTimeDf.laserPower .≈ powerPerGenotype[gen])
    theMean = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => mean)
    time = theMean[:,1]
    ts = (time .< 15) .& (time .>= -5)
    time = time[ts]
    theMean = theMean[ts,2]
    # theMax = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->percentile(p,95)))[ts,2]
    # theMin = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->percentile(p,5)))[ts,2]
    theMax = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->Lensman.plus_sem(p)))[ts,2]
    theMin = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->Lensman.minus_sem(p)))[ts,2]
    meantraceDF = vcat(meantraceDF, DataFrame(mean=theMean, max=theMax, min=theMin,
        time=time, genotype=fill(gen, size(theMean)...)))
end
##

# p = Gadfly.plot(meanTimeDf[idxs,:], x=:time, y=:df_f_mean,
#     color=:cellID, Gadfly.Geom.line, Gadfly.Guide.title("fChRmine"))
meantraceDFnof = filter(x->x.genotype != "fChRmine", meantraceDF)
meantraceDFnof = sort(meantraceDFnof, [:time])
##
# data(meantraceDFnof) * mapping(:time,:mean,col=:genotype) * visual(Lines) |> aog.draw
fontsize_theme = Theme(fontsize = 7 * 300/72, font="Arial")
set_theme!(fontsize_theme)

if "time" in names(meantraceDFnof)
    rename!(meantraceDFnof, :time => "time (s)")
end
meantraceDFnof[meantraceDFnof.genotype .== "ChRmine",:genotype] .= "ChRmine (n=12 cells)"
meantraceDFnof[meantraceDFnof.genotype .== "rsChRmine",:genotype] .= "rsChRmine (n=12 cells)"
d = data(meantraceDFnof)
# p1 = d * mapping(:time,:min,:max,col=:genotype) * visual(Band,alpha=0.3)
W = 183u"mm" - 110u"mm"
H = 50u"mm"
dpi = 300
mm2px = mm->Int(round(uconvert(u"inch", mm) * dpi / 1u"inch" * 0.9))
W = mm2px(W)
H = mm2px(H)
f = Figure(resoultion=(W,H))
p1 = d * mapping("time (s)",:min,:max,col=:genotype) * visual(Band, color="lightblue")
p2 = d * mapping("time (s)",:mean,col=:genotype) * visual(Lines,color="darkblue",linewidth=4)
fig = aog.draw(p1 + p2; axis=(xgridvisible=false, ygridvisible=false,ylabel="Δf/f"))
axs = fig.figure.content[1:2]
map(axs) do ax
    poly!(ax,
    [GeometryBasics.Rect(0,1.8,1,0.1)],
    color = "red",
    colormap = :heat)
end

@show ppath = joinpath(plot_dir,"single-cell-stim")
pt_per_unit = 72/300
save(ppath*".svg", fig, pt_per_unit=pt_per_unit)
save(ppath*".png", fig, pt_per_unit=pt_per_unit)
save(ppath*".pdf", fig, pt_per_unit=pt_per_unit)
fig
