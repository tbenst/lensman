ENV["DISPLAY"] = "localhost:11.0"
##
using DataFrames, Statistics, Arrow, AlgebraOfGraphics, CairoMakie, StatsKit
Data = AlgebraOfGraphics.data
aog = AlgebraOfGraphics

using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Arrow, StatsBase, Unitful, HypothesisTests, GLM
import PyPlot
import Lensman: @pun, @assign

fontsize_theme = Theme(fontsize = 7 * 300/72)
set_theme!(fontsize_theme)
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

major_regions = unique(region_df[!,"major region"])
major_regions = filter(x->x!="Spinal Cord", major_regions)
##
REGION_LIST = [
    # "Tectum Stratum Periventriculare",
    "Optic Tectum",
    "Cerebellum",
    "Habenula",
    "Raphe - Superior",
    # "Raphe - Inferior",
    "Rhombomere 1",
    # "Rhombomere 2",
    # "Rhombomere 3",
    # "Rhombomere 4",
    "Subpallial Gad1b cluster",
    "Gad1b Cluster 1",
    "Gad1b Cluster 2",
    # "Gad1b Cluster 12",
    "Gad1b Cluster 16",
    # "Otpb Cluster 4",
    "Glyt2 Cluster 12",
    "Glyt2 Stripe 3",
]

numcol = 3

# TODO write abstraction for this for any num arguments to pass through
function quick_mean_sem(df_f, region, subregion, period)
    DataFrame("Δf/f"=>mean(df_f), "sem"=>sem(df_f),
    "plus_sem"=>Lensman.plus_sem(df_f), "minus_sem"=>Lensman.minus_sem(df_f),
    "region"=>region,"period"=>period, "subregion"=>subregion)
end

## Second-order proj field mapping all brain regions
second_order_uri = secondorder_names[1]
sdf = filter(row->row[:uri]==second_order_uri, region_df);
plot_df = filter(row -> row.area > 100, sdf)
plot_df = Lensman.rename_regions_in_df(plot_df)

stim_names = sort(unique(plot_df[:,:stim]))
pretty_stim_names = ["habenula\n2nd order", "raphe", "outside\nbrain"]
plot_df[:,"stim target"] = string.(plot_df[:,:stim])
for (m,mm) in zip(stim_names,pretty_stim_names)
    idxs = plot_df[:,:stim] .== m
    plot_df[idxs,"stim target"] .= mm
end

fig = Figure(resolution=(2000, 3000))

for (ii,major) in enumerate(major_regions)
    times, remainder = divrem(ii,numcol)
    i = times + 1
    j = remainder + 1
    pdf = filter(r->r["major region"] == major, plot_df)
    pdf = combine(groupby(pdf, [:region,Symbol("stim target")]),
        [Symbol("Δf/f"), :region, Symbol("stim target")]=>quick_mean_sem=>AsTable)
    # https://discourse.julialang.org/t/makie-errorbars-for-grouped-bar-graphs/62361/7
    # p1 = Data(pdf) * visual(BarPlot) *
    #     mapping(:region, "Δf/f", color="stim target", dodge="stim target")
    p2 = Data(pdf) * visual(CrossBar) *
        mapping(:region, "Δf/f", "minus_sem", "plus_sem", color="stim target", dodge="stim target")
        # mapping(:region, "Δf/f", color="stim target", dodge="stim target")
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
    global grid = aog.draw!(fig[i,j], p, axis=(xticklabelrotation = pi/6,))
    ylims!(-0.2,0.75)
end

legend!(fig[1,end + 1], grid)
ppath = joinpath(plot_dir,
    "all-region_second_order_region_df_f_$(replace(second_order_uri,"""/"""=>"""_"""))")
@show ppath
pt_per_unit = 72/300 # Makie defaults to 72px
save(ppath*".png", fig, pt_per_unit=pt_per_unit)
save(ppath*".svg", fig, pt_per_unit=pt_per_unit)
save(ppath*".pdf", fig, pt_per_unit=pt_per_unit)
fig


## Second-order proj field mapping: select regions

W = 183u"mm" - 110u"mm"
H = 53u"mm"
dpi = 300
mm2px = mm->Int(round(uconvert(u"inch", mm) * dpi / 1u"inch" * 0.9))
W = mm2px(W)
H = mm2px(H)

fig = Figure(resolution=(W, H))
pdf = combine(groupby(plot_df, [:region,:subregion, Symbol("stim target")]),
    [Symbol("Δf/f"), :region, Symbol("stim target")]=>quick_mean_sem=>AsTable)
pdf[pdf.subregion .== "Tectum Stratum Periventriculare",:subregion] .= "Optic Tectum"
pdf = filter(row -> row.region in REGION_LIST, pdf)
for region in REGION_LIST
    @assert sum(pdf.region .== region) > 0 region
end

p2 = Data(pdf) * visual(CrossBar) *
    mapping(:region, "Δf/f", "minus_sem", "plus_sem", color="stim target", dodge="stim target")
grid = aog.draw!(fig, p2, axis=(xticklabelrotation = pi / 6,))
legend!(fig[1,end + 1], grid)
ppath = joinpath(plot_dir,
    "select-region_second_order_region_df_f_$(replace(second_order_uri,"""/"""=>"""_"""))")
@show ppath
pt_per_unit = 72/300 # Makie defaults to 72px
save(ppath*".png", fig, pt_per_unit=pt_per_unit)
save(ppath*".svg", fig, pt_per_unit=pt_per_unit)
save(ppath*".pdf", fig, pt_per_unit=pt_per_unit)
fig

## Titration barplot
# only 3 stimuli...?
# uri = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
uri = "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"
# uri = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
titration_plot_df = filter(row->row.uri == uri,
    region_df);
titration_plot_df = Lensman.add_period_to_df(titration_plot_df);
titration_plot_df = filter(row->row.uri == uri,
    titration_plot_df);
titration_plot_df = Lensman.rename_regions_in_df(titration_plot_df)
# titration_plot_df = filter(row->row.region in mask_names,
#     titration_plot_df);
# fig = Figure(resolution = (12000, 500))
fig = Figure(resolution = (2000, 3000))
# ga = fig[1, 1] = GridLayout()
# TODO: recreate panel b...
["habenula", "optic tectum", "rhombomere 1", "raphe", "cerebellum"]
# ps = []


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


if uri =="2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"
    ymax = 0.1
else
    ymax = 0.3
end

for (ii,major) in enumerate(sort(major_regions))
    times, remainder = divrem(ii,numcol)
    i = times + 1
    j = remainder + 1
    pdf = filter(r->r["major region"] == major, titration_plot_df)
    pdf = combine(groupby(pdf, [:subregion,:period]),
        [Symbol("Δf/f"), :subregion, :period]=>quick_mean_sem=>AsTable)
    # https://discourse.julialang.org/t/makie-errorbars-for-grouped-bar-graphs/62361/7
    # p1 = Data(pdf) * visual(BarPlot) *
    #     mapping(:subregion, "Δf/f", color=:period, dodge=:period)
    p2 = Data(pdf) * visual(CrossBar) *
        mapping(:subregion, "Δf/f", "minus_sem", "plus_sem", color=:period, dodge=:period)
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
    grid = aog.draw!(fig[i,j], p, axis=(xticklabelrotation = pi/6,))
    ylims!(-0.2,0.3)
end

# p1 = Data(titration_plot_df) * visual(BoxPlot, show_notch=true, show_outliers=false) *
# mapping(:region, "Δf/f", color=:period, dodge=:period, row="major region")
# grid = draw!(fig,p1, axis=(xticklabelrotation = pi/6,))
legend!(fig[1,end+1], grid)
ppath = joinpath(plot_dir,
    "all-region_titration_df_f_by-trial_$(replace(uri,"""/"""=>"""_"""))")
@show ppath
save(ppath*".png", fig)
save(ppath*".svg", fig)
save(ppath*".pdf", fig)
fig

## Titration select regions
W = 183u"mm" / 1.5
# W = 183u"mm" / 2
H = 57u"mm"
dpi = 300 # makie default; we scale up on save
mm2px = mm->Int(round(uconvert(u"inch", mm) * dpi / 1u"inch"))
@show W = mm2px(W)
@show H = mm2px(H)
fig = Figure(resolution=(W, H))

pdf = combine(groupby(titration_plot_df, [:region,:period]),
    [Symbol("Δf/f"), :region, :subregion, :period]=>quick_mean_sem=>AsTable)

pdf[pdf.subregion .== "Tectum Stratum Periventriculare",:subregion] .= "Optic Tectum"
pdf = filter(row -> row.subregion in REGION_LIST, pdf)
@warn "we initially chose subregions so no collisions; could erroneously combine subregions if violated"
pdf.period .= replace.(Array(pdf.period), "early" => "1-6")
pdf.period .= replace.(Array(pdf.period), "mid" => "7-11")
pdf.period .= replace.(Array(pdf.period), "late" => "12-16")
rename!(pdf, "period" => "# bilateral\ncell stim")
period = "# bilateral\ncell stim"
p2 = Data(pdf) * visual(CrossBar, width=0.5) *
    mapping(:subregion, "Δf/f", "minus_sem", "plus_sem", color=period, dodge=period)

colors = colormap("Blues",4)[2:4]
grid = aog.draw!(fig, p2, axis=(xticklabelrotation = pi / 5.5,),palettes=(color=colors,))
legend!(fig[1,end + 1], grid)
ppath = joinpath(plot_dir,
    "select-region_titration_df_f_by-trial_$(replace(uri,"""/"""=>"""_"""))")
@show ppath
# pt_per_unit = 72/300 # Makie defaults to 72px
# pt_per_unit = 300/72 # Makie defaults to 72px
pt_per_unit = 72/300 # Makie defaults to 72px
save(ppath*".png", fig, pt_per_unit=pt_per_unit)
save(ppath*".svg", fig, pt_per_unit=pt_per_unit)
save(ppath*".pdf", fig, pt_per_unit=pt_per_unit)
fig

## ANOVA
periods = ["1-6", "7-11", "12-16"]
region_names = unique(titration_plot_df[!,:region])
pvals = map(region_names) do region
    groups = [Float64.(filter(r->(string(r[:period]) == p) && (r.region == region), titration_plot_df)[:,"Δf/f"])
        for p in ["early", "mid", "late"]]
    # KW = KruskalWallisTest(groups...)
    # pvalue(KW)

    df1 = DataFrame(dose = 1, score=groups[1])
    df2 = DataFrame(dose = 2, score=groups[2])
    df3 = DataFrame(dose = 3, score=groups[3])
    glm_df = vcat(df1, df2, df3)


    nullModel = lm(@formula(score ~ 1), glm_df);
    doseModel = lm(@formula(score ~ 1 + dose), glm_df);

    F = ftest(nullModel.model, doseModel.model);
    F.pval[2]
end
bonferoni_p5 = 0.05 / length(pvals)
bonferoni_p1 = 0.01 / length(pvals)
significant5 = pvals .< bonferoni_p5
significant1 = pvals .< bonferoni_p1
println("=====Significant at p=0.05: $(sum(significant5))/$(length(pvals))=====")
for n in sort(region_names[significant5]); println(n); end
println("=====Significant at p=0.01: $(sum(significant1))/$(length(pvals))=====")
for n in sort(region_names[significant1]); println(n); end
# see output at bottom

##

uri

## for 2021-06-08 titration fish
=====Significant at p=0.05: 41/157=====
Diencephalon -
Diencephalon - Anterior pretectum cluster of vmat2 Neurons
Diencephalon - Dorsal Thalamus
Diencephalon - Habenula
Diencephalon - Migrated Area of the Pretectum (M1)
Diencephalon - Olig2 Band
Diencephalon - Olig2 Band 2
Diencephalon - Pineal
Diencephalon - Pineal Vmat2 cluster 
Diencephalon - Pretectal Gad1b Cluster
Diencephalon - Pretectal dopaminergic cluster
Diencephalon - Pretectum
Diencephalon - Ventral Thalamus
Mesencephalon -
Mesencephalon - Oculomotor Nucleus nIII
Mesencephalon - Otpb Cluster
Mesencephalon - Oxtl Cluster Sparse
Mesencephalon - Retinal Arborization Field 8 (AF8)
Mesencephalon - Retinal Arborization Field 9 (AF9)
Mesencephalon - Tecum Neuropil
Mesencephalon - Tegmentum
Mesencephalon - Torus Longitudinalis
Mesencephalon - Vglut2 cluster 1
Rhombencephalon - Gad1b Cluster 1
Rhombencephalon - Gad1b Cluster 17
Rhombencephalon - Gad1b Cluster 2
Rhombencephalon - Isl1 Cluster 1
Rhombencephalon - Isl1 Cluster 2
Rhombencephalon - Lobus caudalis cerebelli
Rhombencephalon - Neuropil Region 6
Rhombencephalon - Oculomotor Nucleus nIV
Rhombencephalon - Otpb Cluster 4
Rhombencephalon - Raphe - Superior
Rhombencephalon - Vglut2 cluster 1
Telencephalon -
Telencephalon - Pallium
Telencephalon - Subpallial Gad1b cluster
Telencephalon - Subpallial dopaminergic cluster
Telencephalon - Subpallium
Telencephalon - Vglut2 rind
Telencephalon - Vmat2 cluster
=====Significant at p=0.01: 38/157=====
Diencephalon -
Diencephalon - Anterior pretectum cluster of vmat2 Neurons
Diencephalon - Dorsal Thalamus
Diencephalon - Habenula
Diencephalon - Migrated Area of the Pretectum (M1)
Diencephalon - Olig2 Band
Diencephalon - Olig2 Band 2
Diencephalon - Pineal
Diencephalon - Pineal Vmat2 cluster 
Diencephalon - Pretectal Gad1b Cluster
Diencephalon - Pretectal dopaminergic cluster
Diencephalon - Pretectum
Diencephalon - Ventral Thalamus
Mesencephalon -
Mesencephalon - Oculomotor Nucleus nIII
Mesencephalon - Otpb Cluster
Mesencephalon - Oxtl Cluster Sparse
Mesencephalon - Retinal Arborization Field 8 (AF8)
Mesencephalon - Retinal Arborization Field 9 (AF9)
Mesencephalon - Tecum Neuropil
Mesencephalon - Tegmentum
Mesencephalon - Torus Longitudinalis
Mesencephalon - Vglut2 cluster 1
Rhombencephalon - Gad1b Cluster 1
Rhombencephalon - Gad1b Cluster 17
Rhombencephalon - Gad1b Cluster 2
Rhombencephalon - Isl1 Cluster 1
Rhombencephalon - Isl1 Cluster 2
Rhombencephalon - Neuropil Region 6
Rhombencephalon - Otpb Cluster 4
Rhombencephalon - Raphe - Superior
Rhombencephalon - Vglut2 cluster 1
Telencephalon -
Telencephalon - Pallium
Telencephalon - Subpallial Gad1b cluster
Telencephalon - Subpallium
Telencephalon - Vglut2 rind
Telencephalon - Vmat2 cluster