ENV["DISPLAY"] = "localhost:11.0"
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

plot_dir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/2021-09-01_regions"
mkdir(plot_dir)

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
##
sort!(df, "Δf/f")
good_idxs = (~).(isnan.(df[:,"Δf/f"]))
df = df[good_idxs,:];
##
# df[:,"one_minus_three"] =

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
sdf = filter(row->row[:uri]==secondorder_names[1], df);

mask_names = [
    "Diencephalon - Habenula",
    "Mesencephalon - Tectum Stratum Periventriculare",  # only one
    "Rhombencephalon - Cerebellum",  # choosing most general one
    "Rhombencephalon - Raphe - Superior",
    "Rhombencephalon - Rhombomere 1",
    "Diencephalon - Habenula",
    "Rhombencephalon - Gad1b Cluster 2"
]
pretty_names = [
    "Habenula",
    "Optic tectum",  # only one
    "Cerebellum",  # choosing most general one
    "Raphe - Superior",
    "Rhombomere 1",
    "Habenula",
    "Gad1b Cluster 2"
]
plot_df = filter(row -> any(mask_names .== row.region), sdf)
unique(plot_df[:, :region])
plot_df[:,:regions] = plot_df[:, :region]
for (m,mm) in zip(mask_names,pretty_names)
    idxs = plot_df[:,:region] .== m
    plot_df[idxs,:regions] .= mm
end
stim_names = sort(unique(plot_df[:,:stim]))
pretty_stim_names = ["habenula second order", "raphe", "outside brain"]
plot_df[:,"stim target"] = plot_df[:,:stim]
for (m,mm) in zip(stim_names,pretty_stim_names)
    idxs = plot_df[:,:stim] .== m
    plot_df[idxs,"stim target"] .= mm
end
##
fig = Figure(resolution=(1200, 500))
p1 = Data(plot_df) * visual(BoxPlot, show_notch=true) *
    mapping(:regions, "Δf/f", color="stim target", dodge="stim target")
grid = draw!(fig, p1, axis=(xticklabelrotation = pi / 6,))
legend!(fig[1,end + 1], grid)
fig
ppath = joinpath(plot_dir, "second_order_region_df_f")
save(ppath*".png", fig)
save(ppath*".svg", fig)
