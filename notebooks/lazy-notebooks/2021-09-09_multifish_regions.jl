ENV["DISPLAY"] = "localhost:12.0"
##
using DataFrames, Statistics, Arrow, AlgebraOfGraphics, CairoMakie
Data = AlgebraOfGraphics.data

using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Arrow
import PyPlot
import Lensman: @pun, @assign
##
# rsync -mrltvP --include "*kalman/regions_df.arrow" --include "*/" --exclude "*" osprey:/scratch/allan/b115/ .

plot_dir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/2021-09-09_regions"
mkpath(plot_dir)

# transferred
arrow_paths = [
    "2021-06-01_rsChRmine_h2b6s/fish3/h5_output/TSeries-IPNraphe-118trial-072_kalman.h5/kalman/regions_df.arrow"
    "2021-06-08_rsChRmine_h2b6s/fish2/h5_output/TSeries-lrhab-titration-123_kalman.h5/kalman/regions_df.arrow"
    "2021-06-08_rsChRmine_h2b6s/fish2/h5_output/TSeries-lrhab-118trial-122_kalman.h5/kalman/regions_df.arrow"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/h5_output/TSeries-lrhab-118trial-061_kalman.h5/kalman/regions_df.arrow"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/h5_output/TSeries-titration-192trial-062_kalman.h5/kalman/regions_df.arrow"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/h5_output/TSeries-lrhab-118trial-069_kalman.h5/kalman/regions_df.arrow"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/h5_output/TSeries-titration-192trial-070_kalman.h5/kalman/regions_df.arrow"
]

recording_names = map(arrow_paths) do a
    p = joinpath(splitpath(a)[[1,2,4]]...)
    replace(p, "_kalman.h5" => "")
end
titration_names = filter(x->occursin("titration",x), recording_names)
secondorder_names = filter(x->occursin("IPN",x), recording_names)
lrhab_names = filter(x->~occursin("titration",x) & occursin("lrhab",x),
    recording_names)

##
root_dir = "/data/dlab/b115"
df = DataFrame()
for (p,r) in zip(arrow_paths, recording_names)
    rdf = Arrow.Table(joinpath(root_dir, p)) |> DataFrame
    rdf[:,:uri] .= r
    df = vcat(df,rdf)
end
if typeof(df.stim[1]) == String
    df[!,:stim] = parse.(Int,df.stim)
end;
df = Lensman.add_major_regions_to_df(df)
# julia> names(df)
# 9-element Vector{String}:
#  "Δf/f"
#  "f"
#  "f0"
#  "area"
#  "region"
#  "stim"
#  "trial"
#  "hemisphere"
#  "uri"

# combine by region
function calc_df_both_hemispheres(f0, area,f,region,df_f,stim,trial,uri,major)
    newf0 = sum(area .* f0)
    newf = sum(area .* f)
    DataFrame("f0"=>newf0, "area"=>sum(area),"f"=>newf, "Δf/f" => (newf-newf0) / newf0,
    "region"=>region[1], "stim"=>stim[1], "trial"=>trial[1], "uri"=>uri[1], "major region"=>major)
end

gb = groupby(df, [:region, :uri, :trial])
@assert size(gb[1],1) == 2 # left&right hemisphere
region_df = combine(gb,
    [:f0, :area, :f,:region,Symbol("Δf/f"),:stim,:trial,:uri, Symbol("major region")]
    => calc_df_both_hemispheres => AsTable)
sort!(df, "Δf/f")
sort!(region_df, "Δf/f")
good_idxs = (~).(isnan.(df[:,"Δf/f"]))
df = df[good_idxs,:];
good_idxs = (~).(isnan.(region_df[:,"Δf/f"]))
region_df = region_df[good_idxs,:];
##
for r in lrhab_names
    rdf = filter(row->(row[:uri]==r) & (row[:stim]=="3") &
        (row[:area] > 300), df)
    # @show rdf[1:20,["Δf/f", "f0", "region"]]
    # @show rdf[end-20:end,["Δf/f", "area"]]
    println((rdf[end-20:end,"region"] .* "\n")...)
    # @show filter(r->r[:region]=="Diencephalon - Habenula",rdf)
end

##
sdf = filter(row->row[:uri]==secondorder_names[1], region_df);

mask_names = [
    "Diencephalon - Habenula",
    "Mesencephalon - Tectum Stratum Periventriculare",  # only one
    "Rhombencephalon - Cerebellum",  # choosing most general one
    # "Rhombencephalon - Raphe - Superior",
    "Rhombencephalon - Rhombomere 1",
    "Diencephalon - Habenula",
    "Rhombencephalon - Gad1b Cluster 1",
    "Rhombencephalon - Gad1b Cluster 2",
    # "Rhombencephalon - Gad1b Cluster 3",
    "Rhombencephalon - Gad1b Cluster 4",
    "Rhombencephalon - Gad1b Cluster 5",
    "Rhombencephalon - Gad1b Cluster 6",
    "Rhombencephalon - Gad1b Cluster 7",
    "Rhombencephalon - Gad1b Cluster 8",
    "Rhombencephalon - Gad1b Cluster 9",
    "Rhombencephalon - Gad1b Cluster 10",
    "Rhombencephalon - Gad1b Cluster 11",
    "Rhombencephalon - Gad1b Cluster 12",
    "Rhombencephalon - Gad1b Cluster 13",
    "Rhombencephalon - Gad1b Cluster 14",
    "Rhombencephalon - Gad1b Cluster 15",
    "Rhombencephalon - Gad1b Cluster 16",
    "Rhombencephalon - Gad1b Cluster 17",
    "Rhombencephalon - Gad1b Cluster 18",
    "Rhombencephalon - Gad1b Cluster 19",
    "Rhombencephalon - Gad1b Cluster 20"
]
pretty_names = [
    "Habenula",
    "Optic tectum",  # only one
    "Cerebellum",  # choosing most general one
    # "Raphe - Superior",
    "Rhombomere 1",
    "Habenula",
    "Gad1b Cluster 1",
    "Gad1b Cluster 2",
    # "Gad1b Cluster 3",
    "Gad1b Cluster 4",
    "Gad1b Cluster 5",
    "Gad1b Cluster 6",
    "Gad1b Cluster 7",
    "Gad1b Cluster 8",
    "Gad1b Cluster 9",
    "Gad1b Cluster 10",
    "Gad1b Cluster 11",
    "Gad1b Cluster 12",
    "Gad1b Cluster 13",
    "Gad1b Cluster 14",
    "Gad1b Cluster 15",
    "Gad1b Cluster 16",
    "Gad1b Cluster 17",
    "Gad1b Cluster 18",
    "Gad1b Cluster 19",
    "Gad1b Cluster 20"
]
@assert length(mask_names) == length(pretty_names)
plot_df = filter(row -> any(mask_names .== row.region), sdf)
plot_df = filter(row -> row.area > 100, plot_df)
unique(plot_df[:, :region])
plot_df[:,:region] = plot_df[:, :region]
for (m,mm) in zip(mask_names,pretty_names)
    idxs = plot_df[:,:region] .== m
    plot_df[idxs,:region] .= mm
end
stim_names = sort(unique(plot_df[:,:stim]))
pretty_stim_names = ["habenula second order", "raphe", "outside brain"]
plot_df[:,"stim target"] = string.(plot_df[:,:stim])
for (m,mm) in zip(stim_names,pretty_stim_names)
    idxs = plot_df[:,:stim] .== m
    plot_df[idxs,"stim target"] .= mm
end

fig = Figure(resolution=(1200, 500))
p1 = Data(plot_df) * visual(BoxPlot, show_notch=true, show_outliers=false)
# p1 *= mapping("stim target", "Δf/f", color="stim target", dodge="stim target")
p1 *= mapping("region", "Δf/f", color="stim target", dodge="stim target")
grid = draw!(fig, p1, axis=(xticklabelrotation = pi / 6,))
legend!(fig[1,end + 1], grid)
ppath = joinpath(plot_dir, "second_order_region_df_f")
@show ppath
fig
# save(ppath*".png", fig)
# save(ppath*".svg", fig)
# save(ppath*".pdf", fig)


## Titration barplot
plot_df = filter(row->row.uri == "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123",
    region_df);
plot_df = Lensman.add_period_to_df(plot_df);
plot_df = filter(row->row.uri == "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123",
    plot_df);
plot_df = Lensman.rename_regions_in_df(plot_df)
# plot_df = filter(row->row.region in mask_names,
#     plot_df);
# fig = Figure(resolution = (12000, 500))
fig = Figure(resolution = (2000, 3000))
# ga = fig[1, 1] = GridLayout()
# TODO: recreate panel b...
["habenula", "optic tectum", "rhombomere 1", "raphe", "cerebellum"]
# ps = []
major_regions = unique(plot_df[!,"major region"])
major_regions = filter(x->x!="Spinal Cord", major_regions)

numcol = 3


# function calc_df_both_hemispheres(f0, area,f,region,df_f,stim,trial,uri,major)
#     newf0 = sum(area .* f0)
#     newf = sum(area .* f)
#     DataFrame("f0"=>newf0, "area"=>sum(area),"f"=>newf, "Δf/f" => (newf-newf0) / newf0,
#     "region"=>region[1], "stim"=>stim[1], "trial"=>trial[1], "uri"=>uri[1], "major region"=>major)
# end

# gb = groupby(df, [:region, :uri, :trial])
# @assert size(gb[1],1) == 2 # left&right hemisphere
# region_df = combine(gb,
#     [:f0, :area, :f,:region,Symbol("Δf/f"),:stim,:trial,:uri, Symbol("major region")]
#     => calc_df_both_hemispheres => AsTable)


function quick_mean_sem(df_f, region, period)
    DataFrame("Δf/f"=>mean(df_f), "sem"=>sem(df_f),
    "plus_sem"=>Lensman.plus_sem(df_f), "minus_sem"=>Lensman.minus_sem(df_f),
    "region"=>region,"period"=>period)
end

for (ii,major) in enumerate(major_regions)
    times, remainder = divrem(ii,numcol)
    i = times + 1
    j = remainder + 1
    pdf = filter(r->r["major region"] == major, plot_df)
    pdf = combine(groupby(pdf, [:region,:period]),
        [Symbol("Δf/f"), :region, :period]=>quick_mean_sem=>AsTable)
    # https://discourse.julialang.org/t/makie-errorbars-for-grouped-bar-graphs/62361/7
    # p1 = Data(pdf) * visual(BarPlot) *
    #     mapping(:region, "Δf/f", color=:period, dodge=:period)
    p2 = Data(pdf) * visual(CrossBar) *
        mapping(:region, "Δf/f", "minus_sem", "plus_sem", color=:period, dodge=:period)
        # mapping(:region, "Δf/f", color=:period, dodge=:period)
    # p = p1 + p2
    p = p2
    # push!(ps, p)
    # add title, but only the title
    Makie.Axis(fig[i,j], title=major, yticksvisible=false, xticksvisible=false,
        rightspinevisible=false,leftspinevisible=false,topspinevisible=false,
        xgridvisible=false,ygridvisible=false,xlabelvisible=false,ylabelvisible=false,
        xminorgridvisible=false,yminorgridvisible=false,
        xminorticksvisible=false,yminorticksvisible=false,
        xticklabelsvisible=false,yticklabelsvisible=false)
    grid = draw!(fig[i,j], p, axis=(xticklabelrotation = pi/6,))
    ylims!(-0.2,0.3)
end

# p1 = Data(plot_df) * visual(BoxPlot, show_notch=true, show_outliers=false) *
# mapping(:region, "Δf/f", color=:period, dodge=:period, row="major region")
# grid = draw!(fig,p1, axis=(xticklabelrotation = pi/6,))
legend!(fig[1,end+1], grid)
ppath = joinpath(plot_dir, "all_region_df_f_by-trial")
@show ppath
save(ppath*".png", fig)
save(ppath*".svg", fig)
save(ppath*".pdf", fig)
fig
