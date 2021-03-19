
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
import Gadfly
using Unitful: μm, m, s

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

# round 2
h33r_path = "$tseriesRootDir/2021-03-16_h33r-chrmine_h2b6s/fish4/TSeries_64cell_8concurrent_2power_8rep-607"
rs_path = "$tseriesRootDir/2021-03-16_rschrmine_h2b6s/fish3/TSeries_64cell_8concurrent_2power_8rep-407"
gcamp_path = "$tseriesRootDir/2021-03-16_h2b6s/fish1/TSeries_64cell_8concurrent_2power_8rep-207"
wt_path = "$tseriesRootDir/2021-03-16_wt-chrmine_h2b6s/fish2/TSeries_64cell_8concurrent_2power_8rep-221"

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
h33r_cells, h33r = readExp(h33r_path, "H33R-ChRmine", 1e4)
rs_cells, rs = readExp(rs_path, "rsChRmine", 2e4)
gcamp_cells, gcamp = readExp(gcamp_path, "control", 2e4)


fluor = vcat(h33r, rs, wt, gcamp);
plotDir = "$tseriesRootDir/../plots"

##

# hack to resample data for same sampling rate...
nTimes = map(x-> length(unique(x.time)), groupby(fluor, :genotype))
maxTime = map(x-> maximum(x.time), groupby(fluor, :genotype))
minTime = map(x-> minimum(x.time), groupby(fluor, :genotype))
frameRates = round.(nTimes ./ (maxTime - minTime), digits=0)
@show frameRates
# resample = Int.(frameRates ./ 3)

# round 1
# resample = [1, 2, 1]
# @warn "resampling hardcode stuff"
# round 2
resample = [1, 1, 1, 1]
# resample = [1, 2,2,2]
# resample = [3, 6,6,6]
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

Gadfly.plot(df_fs, x=:df_f, ygroup=:genotype,
    Gadfly.Geom.subplot_grid(Gadfly.Geom.histogram()))

##
nplots = 32
genotypes = sort(unique(fluorDF.genotype))
@show genotypes
plots = Gadfly.Plot[]
# for gen in genotypes[[1]]

genotype = "ChRmine"
genotype = "rsChRmine"
# genotype = "H33R-ChRmine"
# genotype = "control"
for gen in[genotype]
    cellIDs = unique(fluorDF[fluorDF.genotype .== gen,:cellID])
    nCells = length(cellIDs)
    # idxs = cellIDs[randperm(nCells)[1:nplots]]
    idxs = df_fs[(df_fs.df_f .>= 0.28) .& (df_fs.genotype .== gen), :cellID]
    # idxs = idxs[randperm(length(idxs))[1:nplots]]
    idxs = idxs[1:nplots]
    # @show length(idxs), nCells
    for cellID in idxs
        idxs = meanTimeDf.cellID .== cellID
        # @show sum(idxs)
        push!(plots, Gadfly.plot(meanTimeDf[idxs,:], x=:time, y=:df_f_mean, color=:laserPower,
            Gadfly.Geom.line, Gadfly.Guide.title("cellID $cellID")))
    end
end

# p_nConcurrent = gridstack(reshape(plots[nplots*2+1:nplots*3], 8,5))
p = gridstack(reshape(plots[1:nplots], 8,4))
plotPath = joinpath(plotDir, "$genotype-examples.png")
img = SVG(plotPath, 32inch, 16inch)
Gadfly.draw(img, p)
@show plotPath
p

## round 1
h33r_good_idxs = UInt32[10078, 10087, 10184, 10089, 10216, 10077, 10083, 10124, 10200, 10120, 10192, 10226]
# wt_good_idxs = UInt32[318, 618, 396, 621, 408, 834, 854, 170, 415, 310, 1018, 853]
wt_good_idxs = UInt32[169, 854, 1024, 307, 834, 1005, 310, 1015, 1016, 318, 853, 439] # 51
rs_good_idxs = UInt32[20361, 20341, 20629, 20249, 20641, 20479, 20905, 20960, 20821, 20322, 20373, 20130] #20174


# round 2
# h33r_good_idxs = UInt32[10058, 10005, 10049, 10021, 10047]
control_good_idxs = UInt32[]
wt_good_idxs = UInt32[42, 14, 9, 47, 49, 12, 32, 43, 13, 42, 61, 41, 57 ]
# 20022, 20055 are weak
rs_good_idxs = UInt32[20009, 20010, 20019, 20004, 20032, 20024, 20064, 20020, 20014, 20017, 20022, 20055]

@assert length(unique(h33r_good_idxs)) == 12
@assert length(unique(wt_good_idxs)) == 12
@assert length(unique(rs_good_idxs)) == 12


idxs = vcat(h33r_good_idxs, control_good_idxs, wt_good_idxs, rs_good_idxs)

p
##
toUse = map(x->x in idxs, meanTimeDf.cellID)
line = Gadfly.plot(meanTimeDf[toUse, :], x=:time,
    xgroup=:genotype,
    y="df_f_mean",
    ymin="df_f_minimum",
    ymax="df_f_maximum",
    Gadfly.Geom.subplot_grid(Gadfly.Geom.line, Gadfly.Geom.ribbon))
    # Gadfly.Geom.subplot_grid(Gadfly.Geom.line))
##
line = Gadfly.plot(meanTimeDf[toUse, :], x=:time,
    xgroup=:genotype,
    ygroup=:laserPower,
    y="df_f_mean",
    color=:laserPower,
    Gadfly.Geom.subplot_grid(Gadfly.Geom.line))


##
keys(groupby(fluor, [:genotype, :laserPower]))
keys(groupby(meanTimeDf, [:genotype, :laserPower]))
##

# idxs = map(i->i in h33r_good_idxs, meanTimeDf.cellID)
powerPerGenotype = Dict("ChRmine" => 9.1, "H33R-ChRmine" => 16.0,
    "rsChRmine" => 6.5)
meantraceDF = DataFrame(time=Float64[], mean=Float64[], min=Float64[], max=Float64[], genotype=String[])
for (gen,idxsToUse) in zip(genotypes, [wt_good_idxs, h33r_good_idxs, rs_good_idxs])
    idxs = map(i->i in idxsToUse, meanTimeDf.cellID)
    # idxs = map(i->i in idxsToUse, meanTimeDf.cellID) .&
    #     (meanTimeDf.laserPower .≈ powerPerGenotype[gen])
    theMean = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => mean)
    time = theMean[:,1]
    ts = time .< 15
    time = time[ts]
    theMean = theMean[ts,2]
    theMax = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->percentile(p,95)))[ts,2]
    theMin = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->percentile(p,5)))[ts,2]
    meantraceDF = vcat(meantraceDF, DataFrame(mean=theMean, max=theMax, min=theMin,
        time=time, genotype=fill(gen, size(theMean)...)))
end
##

# p = Gadfly.plot(meanTimeDf[idxs,:], x=:time, y=:df_f_mean,
#     color=:cellID, Gadfly.Geom.line, Gadfly.Guide.title("H33R-ChRmine"))
p = Gadfly.plot(meantraceDF, x=:time, y=:mean, ymin=:min, ymax=:max,
    ygroup=:genotype,
    Gadfly.Geom.subplot_grid(Gadfly.Geom.line, Gadfly.Geom.ribbon, Gadfly.Coord.Cartesian(xmin=-10, xmax=15, ymin=-0.5, ymax=2)),
    # Gadfly.Geom.subplot_grid(),
    Gadfly.Guide.title("ChRmine (12 cells; 95% CI)"),
    Gadfly.Guide.xlabel("time (s)"), Gadfly.Guide.ylabel("Δf/f"))

plotPath = joinpath(plotDir, "all-chrmine.svg")
img = SVG(plotPath, 8inch, 4inch)
Gadfly.draw(img, p)
@show plotPath
p

# TODO: make pretty by resampling time so even; make lines thicker, -10 to 15s