
using Sockets, Observables, Statistics, Images, Lensman,
    Distributions, Unitful, HDF5, Distributed, SharedArrays, Glob,
    CSV, DataFrames, Plots, Dates, ImageDraw, MAT, StatsBase,
    Compose, ImageMagick, Random, PyCall, Arrow, ProgressMeter,
    RollingFunctions
import Gadfly
using Unitful: μm, m, s

## compare ChRmine: WT, RS & H33R
tseriesRootDir = "/oak/stanford/groups/deissero/users/tyler/b115"
# h33r1 = "$tseriesRootDir/2021-02-16_h2b6s_wt-chrmine/fish3/TSeries-256cell-8concurrent-4freq-055"
# h33r1 = h33r1*"_cells_fluorescence.arrow"
# h33r2_path = joinpath(tseriesRootDir, "2021-02-16_6f_h33r_f0_6dpf/fish2/TSeries-256cell-8concurrent-4freq-051")
# h33r2_cells = h33r2_path*"_cellsDF.arrow"
# h33r2 = h33r2_path*"_cells_fluorescence.arrow"
h33r2_path = joinpath(tseriesRootDir, "$tseriesRootDir/2021-02-02_f1_h33r_GC6f_6dpf/fish2/TSeries-1024cell-32concurrent-5power-060")
h33r2_cells = h33r2_path*"_cellsDF.arrow"
h33r2 = h33r2_path*"_cells_fluorescence.arrow"
# rs1_path = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043")
# rs1_cells = rs1_path*"_cellsDF.arrow"
# rs1 = rs1_path*"_cells_fluorescence.arrow"
rs1_path = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish2/TSeries-1024cell-32concurrent-048")
rs1_cells = rs1_path*"_cellsDF.arrow"
rs1 = rs1_path*"_cells_fluorescence.arrow"
# wt1_path = "$tseriesRootDir/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043"
wt1_path = "$tseriesRootDir/2021-02-02_wt_chrmine_GC6f/fish3/TSeries-1024cell-32concurrent-5power-10zplane-077"
wt1_cells = wt1_path*"_cellsDF.arrow"
wt1 = wt1_path*"_cells_fluorescence.arrow"
# tseriesDir = joinpath(tseriesRootDir, "2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2/TSeries-128cell-4concurrent-3power-skip7-044")
h33r2_cells = copy(Arrow.Table(h33r2_cells) |> DataFrame);
rs1_cells = copy(Arrow.Table(rs1_cells) |> DataFrame);
wt1_cells = copy(Arrow.Table(wt1_cells) |> DataFrame);

# h33r1 = copy(Arrow.Table(h33r1) |> DataFrame);
h33r2 = copy(Arrow.Table(h33r2) |> DataFrame);
insertcols!(h33r2, 1+size(h33r2,2), :genotype => "h33r-ChRmine")
insertcols!(h33r2, :stimFreq =>
    [h33r2_cells[findall(h33r2_cells.stimStart .== start)[1],:stimFreq] for start in h33r2.stimStart])
insertcols!(h33r2, :laserPower =>
    [h33r2_cells[findall(h33r2_cells.stimStart .== start)[1],:laserPower] for start in h33r2.stimStart])
h33r2.cellID .+= 10000 # ensure unique

rs1 = copy(Arrow.Table(rs1) |> DataFrame);
insertcols!(rs1, 1+size(rs1,2), :genotype => "rsChRmine")
insertcols!(rs1, :stimFreq =>
    [rs1_cells[findall(rs1_cells.stimStart .== start)[1],:stimFreq] for start in rs1.stimStart])
insertcols!(rs1, :laserPower =>
    [rs1_cells[findall(rs1_cells.stimStart .== start)[1],:laserPower] for start in rs1.stimStart])
rs1.cellID .+= 20000

wt1 = copy(Arrow.Table(wt1) |> DataFrame);
insertcols!(wt1, 1+size(wt1,2), :genotype => "ChRmine")
# insertcols!(wt1, :stimFreq =>
#     [wt1_cells[findall(wt1_cells.stimStart .== start)[1],:stimFreq] for start in wt1.stimStart])
# insertcols!(wt1, :laserPower =>
#     [wt1_cells[findall(wt1_cells.stimStart .== start)[1],:laserPower] for start in wt1.stimStart])


fluor = vcat(h33r2, rs1, wt1);
plotDir = "$tseriesRootDir/../plots"

##

function rawToDf_f(fluor::DataFrame)
    fluorDF = []
    # operate over each cell's trace...
    @warn "hack hardcoded `before`"
    for df in groupby(fluor, [:cellID, :stimStart])
        if df.genotype[1] == "h33r-ChRmine"
            before = 60 # 2 seconds
        else
            before = 6
        end
        kalmanFilt = imageJkalmanFilter(medianfilt(df.f,3))
        f0 = mean(kalmanFilt[1:before])
        df = copy(df)
        df_f = @. (kalmanFilt - f0) / (f0 + 10)
        insertcols!(df, size(df,2), :df_f => df_f)
        select!(df, Not(:f))
        push!(fluorDF, df)
    end
    fluorDF = vcat(fluorDF...)
end

fluorDF = rawToDf_f(fluor)
meanTimeDf = combine(groupby(fluorDF, [:cellID, :time, :stimFreq, :laserPower, :genotype]), :df_f => mean)

##
df_fs = DataFrame(cellID=UInt32[], stimFreq=Int64[], laserPower=Float64[],
    genotype=String[], df_f=Float64[])
for df in groupby(meanTimeDf, [:laserPower, :stimFreq, :cellID])
    cellID, stimFreq, laserPower, genotype = df[1,[:cellID, :stimFreq, :laserPower, :genotype]]
    idxs = (df.time .> 1.1) .& (df.time .<= 3.1)
    push!(df_fs, (cellID=cellID, stimFreq=stimFreq, laserPower=laserPower,
        genotype=genotype, df_f=mean(df.df_f_mean[idxs])))
end

Gadfly.plot(df_fs, x=:df_f, ygroup=:genotype,
    Gadfly.Geom.subplot_grid(Gadfly.Geom.histogram()))

##
nplots = 40
genotypes = unique(fluorDF.genotype)
@show genotypes
plots = Gadfly.Plot[]
for gen in unique(fluor.genotype[[3]])
    cellIDs = unique(fluorDF[fluorDF.genotype .== gen,:cellID])
    nCells = length(cellIDs)
    # idxs = cellIDs[randperm(nCells)[1:nplots]]
    # idxs = df_fs[(df_fs.df_f .>= 0.50) .& (df_fs.genotype .== gen), :cellID]
    idxs = df_fs[(df_fs.df_f .>= 0.75) .& (df_fs.genotype .== gen), :cellID]
    idxs = idxs[randperm(length(idxs))[1:nplots]]
    for cellID in idxs
        idxs = meanTimeDf.cellID .== cellID
        push!(plots, Gadfly.plot(meanTimeDf[idxs,:], x=:time, y=:df_f_mean, color=:laserPower,
            Gadfly.Geom.line, Gadfly.Guide.title("cellID $cellID")))
    end
end

p_nConcurrent = gridstack(reshape(plots[nplots*2+1:nplots*3], 8,5))
plotPath = joinpath(plotDir, "debug.png")
img = SVG(plotPath, 32inch, 16inch)
Gadfly.draw(img, p_nConcurrent)
@show plotPath
p_nConcurrent


h33r_good_idxs = UInt32[10078, 10087, 10184, 10089, 10216, 10077, 10083, 10124, 10200, 10120, 10192, 10226]
@assert length(unique(h33r_good_idxs)) == 12
##

idxs = map(i->i in h33r_good_idxs, meanTimeDf.cellID)
theMean = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => mean)
time = theMean[:,1]
ts = time .< 15
time = time[ts]
theMean = theMean[ts,2]
theMax = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->percentile(p,95)))[ts,2]
theMin = combine(groupby(meanTimeDf[idxs,:], :time), :df_f_mean => (p->percentile(p,5)))[ts,2]
# p = Gadfly.plot(meanTimeDf[idxs,:], x=:time, y=:df_f_mean,
#     color=:cellID, Gadfly.Geom.line, Gadfly.Guide.title("H33R-ChRmine"))
p = Gadfly.plot(x=time, y=theMean, ymin=theMin, ymax=theMax,
    Gadfly.Geom.line, Gadfly.Geom.ribbon, Gadfly.Guide.title("H33R-ChRmine (12 cells; 95% CI)"),
    Gadfly.Guide.xlabel("time (s)"), Gadfly.Guide.ylabel("Δf/f"))

plotPath = joinpath(plotDir, "H33R.svg")
img = SVG(plotPath, 8inch, 4inch)
Gadfly.draw(img, p)
@show plotPath
p