
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
h33r2_path = joinpath(tseriesRootDir, "2021-02-16_6f_h33r_f0_6dpf/fish2/TSeries-256cell-8concurrent-4freq-051")
h33r2_cells = h33r2_path*"_cellsDF.arrow"
h33r2 = h33r2_path*"_cells_fluorescence.arrow"
rs1_path = joinpath(tseriesRootDir, "2021-01-26_rsChRmine_6f_7dpf/fish1/TSeries-31concurrent-168trial-3rep-4power-043")
rs1_cells = rs1_path*"_cellsDF.arrow"
rs1 = rs1_path*"_cells_fluorescence.arrow"
wt1_path = "$tseriesRootDir/2021-01-19_chrmine_kv2.1_6f_7dpf/fish1_chrmine/TSeries-1024cell-32concurrent-4power-043"
wt1_cells = wt1_path*"_cellsDF.arrow"
wt1 = wt1_path*"_cells_fluorescence.arrow"
# tseriesDir = joinpath(tseriesRootDir, "2021-02-23_rsChRmine_f0_h2b6s_6dpf/fish2/TSeries-128cell-4concurrent-3power-skip7-044")
h33r2_cells = copy(Arrow.Table(h33r2_cells) |> DataFrame);
rs1_cells = copy(Arrow.Table(rs1_cells) |> DataFrame);
wt1_cells = copy(Arrow.Table(wt1_cells) |> DataFrame);

# h33r1 = copy(Arrow.Table(h33r1) |> DataFrame);
h33r2 = copy(Arrow.Table(h33r2) |> DataFrame);
insertcols!(h33r2, 1+size(h33r2,2), :genotype => "h33r-ChRmine")
h33r2.cellID .+= 10000 # ensure unique
rs1 = copy(Arrow.Table(rs1) |> DataFrame);
insertcols!(rs1, 1+size(rs1,2), :genotype => "rsChRmine")
rs1.cellID .+= 20000
wt1 = copy(Arrow.Table(wt1) |> DataFrame);
insertcols!(wt1, 1+size(wt1,2), :genotype => "ChRmine")




fluor = vcat(h33r2, rs1, wt1);

## plot next-time-step df histogram
diffs = Float64[]
for df in groupby(fluor, [:cellID, :stimStart])
    diffs = vcat(diffs,diff(df.f))
end
histogram(diffs, yaxis=:log)
## cellID is not matched in two experiments :/
# need to match by (x,y)
xyToCellID = Dict()
for cell in eachrow(cells32)
    x, y, cellID = cell[[:x,:y,:cellID]]
    xyToCellID[(x,y)] = cellID
end
length(xyToCellID)

cellIDmap = Dict()
for i in 1:size(cells8,1)
    x, y, cellID = cells8[i,[:x,:y,:cellID]]
    newID = xyToCellID[(x,y)]
    cellIDmap[cellID] = newID
    cells8[i,:cellID] = newID
end

for i in 1:size(fluor8,1)
    fluor8[i,:cellID] = cellIDmap[fluor8[i,:cellID]]
end

## add stim freq to fluor...

insertcols!(fluor8, :stimFreq =>
    [cells8[findall(cells8.stimStart .== start)[1],:stimFreq] for start in fluor8.stimStart])
insertcols!(fluor32, :stimFreq =>
    [cells32[findall(cells32.stimStart .== start)[1],:stimFreq] for start in fluor32.stimStart])
##

function medianfilter(trace, window=3)
    mapwindow(median, trace, (3,))
end

function rawToDf_f(fluor::DataFrame, before=before)
    fluorDF = []
    # operate over each cell's trace...
    for df in groupby(fluor, [:cellID, :stimStart])
        kalmanFilt = imageJkalmanFilter(rollmedian(df.f,3))
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
insertcols!(fluorDF, :nConcurrent => 8)
meanTimeDf = combine(groupby(fluorDF, [:cellID, :time, :nConcurrent]), :df_f => mean)


## 32 concurrent vs 8 concurrent
# idxs = map(c->c in cellIDrankings[1:nplots], meanTimeDf[meanTimeDf.cellID)
plots = Gadfly.Plot[]
cellIDs = unique(cells8.cellID)
nCells = length(cellIDs)
for cellID in cellIDs[randperm(nCells)[1:nplots]]
    idxs = meanTimeDf.cellID .== cellID
    push!(plots, Gadfly.plot(meanTimeDf[idxs,:], x=:time, y=:df_f_mean, color=:nConcurrent,
        Gadfly.Geom.line, ))
end

p_nConcurrent = gridstack(reshape(plots, 8,5))
plotPath = joinpath(plotDir, "$(expName)_nConcurrent_shuffle_3.png")
img = PNG(plotPath, 32inch, 16inch)
Gadfly.draw(img, p_nConcurrent)
@show plotPath
p_nConcurrent

##
plots = Gadfly.Plot[]
for cellID in cellIDs[randperm(nCells)[1:nplots]]
    idxs = fluorDF.cellID .== cellID
    push!(plots, Gadfly.plot(fluorDF[idxs,:], x=:time, y=:df_f, color=:stimFreq,
        Gadfly.Geom.line, ))
end

p_nConcurrent = gridstack(reshape(plots, 8,5))
plotPath = joinpath(plotDir, "$(expName)_stimFreq.png")
img = PNG(plotPath, 32inch, 16inch)
Gadfly.draw(img, p_nConcurrent)
@show plotPath
p_nConcurrent

## TODO: fit/calculate tau off & plot across population...