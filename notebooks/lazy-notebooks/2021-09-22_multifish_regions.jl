# mostly, we use the previously saved .arrow files from osprey
# but, we also supplement with dorsal/ventral habenula as defined by
# center of mass
# ENV["DISPLAY"] = "localhost:11.0"
##
using DataFrames, Statistics, Arrow, AlgebraOfGraphics, CairoMakie, StatsKit, Makie
Data = AlgebraOfGraphics.data
aog = AlgebraOfGraphics
using ImageView

using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Statistics, PyCall, ProgressMeter, HDF5, Distributed,
      Arrow, StatsBase, Unitful, HypothesisTests, GLM, MultipleTesting
import PyPlot
import Lensman: @pun, @assign
L = Lensman
fontsize_theme = Theme(fontsize = 7 * 300/72)
set_theme!(fontsize_theme)
init_workers()
## we will add the dorsal / ventral habenula by splitting center of mass
resources = Resources();
##
recording_2021_06_08_fish2_titration = Recordings[
    "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
](;resources...);
recording_2021_06_02_fish3_ipn = Recordings[
    "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
](;resources...);
##

function dorsal_ventral_hab_df(recording)
    @pun (region_masks_h5, imaging2zseries_plane, zbrain_mask_names,
        tseries, window_len,
        stim_start_idx, stim_end_idx,
        trial_order, tseriesH, tseriesW, tseriesZ

    ) = recording;

    habenula_name, habenula = L.read_first_mask(region_masks_h5, zbrain_mask_names,
        imaging2zseries_plane,"Habenula")
    lhemisphere_name, lhemisphere = L.read_first_mask(region_masks_h5, zbrain_mask_names,
        imaging2zseries_plane,"Left hemisphere")
    rhemisphere_name, rhemisphere = L.read_first_mask(region_masks_h5, zbrain_mask_names,
        imaging2zseries_plane,"Right hemisphere")

    hab = read_registered_mask(region_masks_h5, habenula_name);
    left_hemi = read_registered_mask(region_masks_h5, "Left hemisphere");
    right_hemi = read_registered_mask(region_masks_h5, "Right hemisphere");

        
    prof = sum(hab,dims=[1,2])[1,1,:]
    cumprof = cumsum(prof)
    halfway = searchsortedfirst(cumprof, cumprof[end] / 2)


    # ymax = maximum(prof)
    # df2 = (z=[halfway, halfway], area=[0, ymax])
    # fig = Figure()
    # ax = Makie.Axis(fig[1, 1], title="dorsal/ventral hab cutoff")
    # p1 = (aog.data((z=1:length(prof), area=prof)) + aog.data(df2) ) * mapping(:z, :area)* visual(Lines)
    # aog.draw!(fig[1,1], p1)
    # fig
    dorsal_hab = copy(hab)
    ventral_hab = copy(hab)
    dorsal_hab[:,:,halfway+1:end] .= false
    ventral_hab[:,:,1:halfway] .= false
    # im = RGB.(dorsal_hab)
    # channelview(im)[[2,3],:,:,:] .= 0
    # channelview(im)[2,:,:,:] = ventral_hab
    # imshow(im)
    masks = [dorsal_hab .& left_hemi, dorsal_hab .& right_hemi,
        ventral_hab .& left_hemi, ventral_hab .& right_hemi]
    masks = map(m->imresize(m, tseriesH, tseriesW, tseriesZ) .> 0, masks);
    region_names = ["left Diencephalon - Dorsal Habenula", "right Diencephalon - Dorsal Habenula",
        "left Diencephalon - Ventral Habenula", "right Diencephalon - Ventral Habenula"]
    regions_df = L.per_trial_regions_df(
        tseries, window_len,
        stim_start_idx, stim_end_idx,
        # ss, se,
        trial_order, masks, region_names);
    return regions_df
end
dorsal_ventral_hab_df_2021_06_08_fish2_titration = dorsal_ventral_hab_df(recording_2021_06_08_fish2_titration)
dorsal_ventral_hab_df_2021_06_08_fish2_titration[:,:uri] .= "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
dorsal_ventral_hab_df_2021_06_02_fish3_ipn = dorsal_ventral_hab_df(recording_2021_06_02_fish3_ipn)
dorsal_ventral_hab_df_2021_06_02_fish3_ipn[:,:uri] .= "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
@assert length(unique(dorsal_ventral_hab_df_2021_06_08_fish2_titration.stim)) == 16
@assert length(unique(dorsal_ventral_hab_df_2021_06_02_fish3_ipn.stim)) == 3
##
# rsync -mrltvP --include "*kalman/regions_df.arrow" --include "*/" --exclude "*" osprey:/scratch/allan/b115/ .

plot_dir = "/home/tyler/Dropbox/Science/manuscripts/2021_chrmine-structure/2021-09-09_regions"
mkpath(plot_dir)

# transferred
arrow_paths = [
    # "2021-06-01_rsChRmine_h2b6s/fish3/h5_output/TSeries-IPNraphe-118trial-072_kalman.h5/kalman/regions_df.arrow" # old df, no exclusion of targets
    "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072_kalman/regions_df.arrow" # only df here with exclusion of targets
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


## compare the two arrows for different exclusion criteria...
# TODO: check if area different before & after / if looks okay...?
DataFrames.combine(groupby(filter(r->(r.uri == "2021-06-02_rsChRmine-h2b6s/fish2/regions_df.arrow") && (r.fullname == "Rhombencephalon - Raphe - Superior"), df), [:stim, :hemisphere]), :area => mean)

## sort regions from rostral to caudal
@pun (region_masks_h5, zbrain_mask_names, zbrain_masks, 
) = recording_2021_06_02_fish3_ipn;

rhom_names = ["Rhombencephalon - Rhombomere $i" for i in 1:7]
rhom_masks = [L.read_mask(zbrain_masks, n)
    for n in rhom_names];

##
function is_subregion(mask1, mask2; tol=0.9)
    sum(mask1 .& mask2) / sum(mask1) > tol
end

function get_major_subregion(name, mask)
    if name == "Spinal Cord"
        major_region = "Spinal Cord"
        subregion = "Spinal Cord"
    else
        m = match(r"^([\w\s]+) -\s?(.*)", name)
        if ~isnothing(m)
            major_region, subregion = m[1], m[2]
        else
            error("no match for $(name)")
        end
    end
    major_region, subregion
end

"missing dorsal / ventral hab...?"
function region_to_major_sub(zbrain_mask_names, zbrain_masks, rhom_names, rhom_masks)
    # since so many more regions in hindbrain, we seperate out major regions...
    region_xloc = DataFrame()
    @showprogress for region in zbrain_mask_names
        mask = read_mask(zbrain_masks,region)
        xs = sum(Array(mask),dims=(1,3))[1,:,1]
        # calculate center of mass
        tot_mass = sum(xs)
        x = sum(collect(1:length(xs)) .* xs) / tot_mass

        # find major region
        in_rhombomere = map(m2->is_subregion(mask,m2), rhom_masks)
        if occursin("Rhombencephalon", region) && any(in_rhombomere)
            major = rhom_names[findfirst(in_rhombomere)]
            _, subregion = get_major_subregion(region, mask)
        else
            major, subregion = get_major_subregion(region, mask)
        end

        push!(region_xloc, (fullname=region, region=major, subregion=subregion, x=x))
    end
    idx = findfirst(region_xloc.fullname .== "Diencephalon - Habenula")
    x,region = region_xloc[idx,[:x,:region]]
    push!(region_xloc, (fullname="Diencephalon - Dorsal Habenula", region=region,
        subregion="Dorsal Habenula", x=x))
    push!(region_xloc, (fullname="Diencephalon - Ventral Habenula", region=region,
        subregion="Ventral Habenula", x=x))
    sort!(region_xloc,:x)

end

# skip left / right hemisphere
regions = region_to_major_sub(zbrain_mask_names[1:end-2], zbrain_masks, rhom_names, rhom_masks)
region_levels = [
    "Ganglia",
    "Spinal Cord",
    "Rhombencephalon - Rhombomere 7",
    "Rhombencephalon - Rhombomere 6",
    "Rhombencephalon - Rhombomere 5",
    "Rhombencephalon - Rhombomere 4",
    "Rhombencephalon - Rhombomere 3",
    "Rhombencephalon - Rhombomere 2",
    "Rhombencephalon - Rhombomere 1",
    "Rhombencephalon",
    "Mesencephalon",
    "Diencephalon",
    "Telencephalon"
]
regions.region = CategoricalArray(regions.region; levels=region_levels)
[(reg, sum(regions.region .== reg)) for reg in levels(regions.region)]
##
original_df = DataFrame()
root_dir = "/data/dlab/b115"
for (p,r) in zip(arrow_paths, recording_names)
    rdf = Arrow.Table(joinpath(root_dir, p)) |> DataFrame
    rdf[:,:uri] .= replace(r, "2021-06-01_rsChRmine_h2b6s/fish3" => "2021-06-02_rsChRmine-h2b6s/fish2")
    original_df = vcat(original_df,rdf)
end
##
if sum(original_df.region .== "Diencephalon - Dorsal Habenula") == 0
    print
    df = vcat(original_df, dorsal_ventral_hab_df_2021_06_08_fish2_titration,
        dorsal_ventral_hab_df_2021_06_02_fish3_ipn)
else
    df = original_df
end

function is_not_bad_region(name)
    ~((name =="right Diencephalon - Dorsal Habenula") || (name =="left Diencephalon - Dorsal Habenula") ||
    (name =="right Diencephalon - Ventral Habenula") || (name =="left Diencephalon - Ventral Habenula"))
end
df = filter(r->is_not_bad_region(r.region), df)

if typeof(df.stim[1]) == String
    df[!,:stim] = parse.(Int,df.stim)
end;
DataFrames.rename!(df, :region => :fullname)

# df = Lensman.add_major_regions_to_df(df)
# df = Lensman.rename_regions_in_df(df)
# okay (no subregion yet)
##
function fullname_to_major_subregion(fullname, regions)
    # no idea why these lines of code exist...? seem to break things
    # # support our new regions...
    # if (fullname == "Diencephalon - Dorsal Habenula") || 
    #    (fullname == "Diencephalon - Ventral Habenula")
    #    fullname = "Diencephalon - Habenula"
    # end
    idx = findfirst(regions.fullname .== fullname)
    if isnothing(idx)
        println("FULLNAME NOT FOUND: $fullname")
    end
    (region=regions[idx,:region], subregion=regions[idx,:subregion])
    # regions[idx,[:region, :subregion]]
end

_ftms = f -> fullname_to_major_subregion(f, regions)
ftms  = f -> _ftms.(f)
DataFrames.transform!(df, :fullname => ftms => AsTable)
##

# set levels so we get desired plotting order
df_region = CategoricalArray(df.region, ordered=true)
levels!(df_region, region_levels)
df[!,:region] = df_region
fullnames =  CategoricalArray(df.fullname; ordered=true)
levels!(fullnames, regions.fullname)
df.fullname = fullnames

df
# okay (dorsal hab still here)
##
df[df.subregion .== "Tectum Stratum Periventriculare",:subregion] .= "Optic Tectum"
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

function calc_df_both_hemispheres(f0,area,f,df_f)
    newf0 = sum(area .* f0)
    newf = sum(area .* f)
    DataFrame("f0"=>newf0, "area"=>sum(area),"f"=>newf, "Δf/f" => (newf-newf0) / newf0)
end

gb = groupby(df, [:fullname, :stim, :trial, :uri, :region, :subregion])
# todo: check for a known region where should have two hemisphere..?
# @assert size(gb[1],1) == 2 # left&right hemisphere 
region_df = DataFrames.combine(gb,
    [:f0, :area, :f,Symbol("Δf/f")] => calc_df_both_hemispheres => AsTable;
    keepkeys=true)
sort!(df, "Δf/f")
sort!(region_df, "Δf/f")
good_idxs = (~).(isnan.(df[:,"Δf/f"]))
df = df[good_idxs,:];
good_idxs = (~).(isnan.(region_df[:,"Δf/f"]))
region_df = region_df[good_idxs,:];
major_regions = levels(region_df[!,"region"])
subregions = levels(region_df[!,"subregion"]) # TODO: should sort in loop for each region
# major_regions = filter(x->x!="Spinal Cord", major_regions)
##
REGION_LIST = [
    "Dorsal Habenula",
    "Ventral Habenula",
    "Hindbrain Glyt2 Stripe 3",
    "Raphe - Superior",
    "Hindbrain Gad1b Cluster 16",
    "Hindbrain Gad1b Cluster 2",
    # "Tectum Stratum Periventriculare",
    # "Optic Tectum",
    # "Cerebellum",
    # "Habenula",
    # "Raphe - Inferior",
    # "Rhombomere 1",
    # "Rhombomere 2",
    # "Rhombomere 3",
    # "Rhombomere 4",
    # "Subpallial Gad1b cluster",
    # "Gad1b Cluster 1",
    # "Gad1b Cluster 12",
    # "Otpb Cluster 4",
    # "Glyt2 Cluster 12",
]

numcol = 3

# TODO write abstraction for this for any num arguments to pass through
function quick_mean_sem(df_f; min=0.02)
    thesem = sem(df_f)
    themean = mean(df_f)
    minsem = maximum([thesem, min])
    DataFrame("Δf/f"=>mean(df_f), "sem"=>thesem,
    "plus_sem"=>themean+minsem, "minus_sem"=>themean-minsem)
end
quick_mean_sem(ones(5))
# function quick_mean_sem(df_f, region, subregion, period)
#     DataFrame("Δf/f"=>mean(df_f), "sem"=>sem(df_f),
#     "plus_sem"=>Lensman.plus_sem(df_f), "minus_sem"=>Lensman.minus_sem(df_f),
#     "region"=>region,"period"=>period, "subregion"=>subregion)
# end

####################################################
## Second-order proj field mapping all brain regions
####################################################
# second_order_uri = secondorder_names[1]

# second_order_uri = "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-IPNraphe-118trial-072"
second_order_uri = "2021-06-02_rsChRmine-h2b6s/fish2/regions_df.arrow" # TODO: fix the URI
sdf = filter(row->row[:uri]==second_order_uri, region_df);
second_order_plot_df = filter(row -> row.area > 100, sdf)
# second_order_plot_df = Lensman.rename_regions_in_df(second_order_plot_df)

stim_names = sort(unique(second_order_plot_df[:,:stim]))
pretty_stim_names = ["habenula\n1st order", "raphe", "outside\nbrain"]
second_order_plot_df[!,"stim target"] = string.(second_order_plot_df[:,:stim])

for (m,mm) in zip(stim_names,pretty_stim_names)
    idxs = second_order_plot_df[:,:stim] .== m
    second_order_plot_df[idxs,"stim target"] .= mm
end

stim_targets = CategoricalArray(second_order_plot_df[!,"stim target"], ordered=true)
levels!(stim_targets, ["outside\nbrain", "raphe","habenula\n1st order"])
second_order_plot_df[!,"stim target"] = stim_targets


ncol_per_region = DataFrames.combine(groupby(second_order_plot_df, :region), :subregion => x->length(unique(x)))
# let's group:
# Rhombomere 7, Spinal Cord, Ganglia
# Rhombomere 1 - 6
# Rhombencephalon
# Mesencephalon, Diencephalon, Telencephalon
# 
##
region_groups = [
    ["Ganglia", "Spinal Cord", "Rhombencephalon - Rhombomere 7"],
    ["Rhombencephalon - Rhombomere $i" for i in 6:-1:1],
    ["Rhombencephalon"],
    ["Mesencephalon", "Diencephalon", "Telencephalon"]
]
region_group_names = replace.(join.(region_groups,", "), "Rhombencephalon - " => "")
W = 183u"mm"
H = 247u"mm"
dpi = 300 # makie default; we scale up on save
mm2px = mm->Int(round(uconvert(u"inch", mm) * dpi / 1u"inch"))
@show W = mm2px(W)
@show H = mm2px(H)

tiny_fontsize_theme = Theme(fontsize = 5 * 300/72)

with_theme(tiny_fontsize_theme) do

    fig = Figure(resolution=(W, H))
    for (ii,(reg_name,reg_group)) in enumerate(zip(region_group_names, region_groups))
        # times, remainder = divrem(ii,numcol)
        # i = times + 1
        # j = remainder + 1
        pdf = filter(r->r[:region] in reg_group, second_order_plot_df)
        pdf = DataFrames.combine(groupby(pdf, [:subregion,Symbol("stim target")]),
            Symbol("Δf/f")=>quick_mean_sem=>AsTable; keepkeys=true)
        # https://discourse.julialang.org/t/makie-errorbars-for-grouped-bar-graphs/62361/7
        # p1 = Data(pdf) * visual(BarPlot) *
        #     mapping(:region, "Δf/f", color="stim target", dodge="stim target")
        p2 = Data(pdf) * visual(CrossBar) *
            mapping(:subregion, "Δf/f", "minus_sem", "plus_sem", color="stim target", dodge="stim target")
            # mapping(:region, "Δf/f", color="stim target", dodge="stim target")
        # p = p1 + p2
        p = p2
        # push!(ps, p)
        # add title, but only the title
        # Makie.Axis(fig[ii,2:9], yticksvisible=false, xticksvisible=false,
        Makie.Axis(fig[ii,1], title=reg_name, yticksvisible=false, xticksvisible=false,
            rightspinevisible=false,leftspinevisible=false,topspinevisible=false,
            xgridvisible=false,ygridvisible=false,xlabelvisible=false,ylabelvisible=false,
            xminorgridvisible=false,yminorgridvisible=false,
            xminorticksvisible=false,yminorticksvisible=false,
            xticklabelsvisible=false,yticklabelsvisible=false)
        global grid = aog.draw!(fig[ii,1], p, axis=(xticklabelrotation = pi/6,))
        ylims!(-0.2,1.5)
    end

    # legend!(fig[1,end + 1], grid)
    ppath = joinpath(plot_dir,
        "all-region_second_order_region_df_f_$(replace(second_order_uri,"""/"""=>"""_"""))")
    @show ppath
    pt_per_unit = 72/300 # Makie defaults to 72px
    # save(ppath*".png", fig, pt_per_unit=pt_per_unit)
    # save(ppath*".svg", fig, pt_per_unit=pt_per_unit)
    # save(ppath*".pdf", fig, pt_per_unit=pt_per_unit)
    fig
end
##################################################
## Second-order proj field mapping: select regions
##################################################
# W = 183u"mm" - 90u"mm"
W = 75u"mm"
# H = 53u"mm"
H = 80u"mm"
dpi = 300
mm2px = mm->Int(round(uconvert(u"inch", mm) * dpi / 1u"inch"))
W = mm2px(W)
H = mm2px(H)

fig = Figure(resolution=(W, H))
pdf = DataFrames.combine(groupby(second_order_plot_df, [:region,:subregion, Symbol("stim target")]),
    Symbol("Δf/f")=>quick_mean_sem=>AsTable, keepkeys=true)
pdf[pdf.subregion .== "Tectum Stratum Periventriculare",:subregion] .= "Optic Tectum"
pdf[pdf.subregion .== "Glyt2 Stripe 3",:subregion] .= "Hindbrain Glyt2 Stripe 3"
pdf[pdf.subregion .== "Gad1b Cluster 16",:subregion] .= "Hindbrain Gad1b Cluster 16"
pdf[pdf.subregion .== "Gad1b Cluster 2",:subregion] .= "Hindbrain Gad1b Cluster 2"

# pdf[pdf[:,"stim target"] .== "outside\nbrain","stim target"] .= "control"
pdf = filter(row -> row.subregion in REGION_LIST, pdf)
# sort!(pdf, "stim target"; lt=(a,b)->CategoricalValue(a,stim_targets)<CategoricalValue(b,stim_targets))
pdf.subregion = CategoricalArray(pdf.subregion, ordered=true)
levels!(pdf.subregion, REGION_LIST)

stim_targets = Array(pdf[:,"stim target"])
stim_targets .= replace.(stim_targets, "\n" => " ")
# stim_targets = CategoricalArray(pdf[:,"stim target"], ordered=true)
stim_targets = CategoricalArray(stim_targets, ordered=true)

levels!(stim_targets, ["outside brain", "raphe","habenula 1st order"])
# levels!(stim_targets, ["outside\nbrain", "raphe","habenula\n1st order"])
pdf[!,"stim target"] = stim_targets

for region in REGION_LIST
    @assert sum(pdf.subregion .== CategoricalValue(region, pdf.subregion)) > 0 region
end

p2 = Data(pdf) * visual(CrossBar, width=0.5) *
    mapping(:subregion, "Δf/f", "minus_sem", "plus_sem", color="stim target", dodge="stim target")
grid = aog.draw!(fig[1:3,2:20], p2, axis=(xticklabelrotation = pi / 6,))
legend!(fig[4,1:20], grid, orientation = :horizontal)
ylims!(-0.15,1.7)
ppath = joinpath(plot_dir,
    "select-region_second_order_region_df_f_$(replace(second_order_uri,"""/"""=>"""_"""))")
@show ppath*".pdf"
pt_per_unit = 72/300 # Makie defaults to 72px
save(ppath*".png", fig, pt_per_unit=pt_per_unit)
save(ppath*".svg", fig, pt_per_unit=pt_per_unit)
save(ppath*".pdf", fig, pt_per_unit=pt_per_unit)
fig

####################################
## Titration barplot
####################################
# only 3 stimuli...?
# uri = "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
# uri = "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-titration-192trial-070"
uri = "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
titration_plot_df = filter(row->row.uri == uri,
    region_df);

titration_plot_df = Lensman.add_period_to_df(titration_plot_df);
titration_plot_df = filter(row->row.uri == uri,
    titration_plot_df);
# titration_plot_df = Lensman.rename_regions_in_df(titration_plot_df)

_ftms = f -> fullname_to_major_subregion(f, regions)
ftms  = f -> _ftms.(f)
DataFrames.transform!(titration_plot_df, :fullname => ftms => AsTable)
# titration_plot_df = filter(row->row.region in mask_names,
#     titration_plot_df);
fig = Figure(resolution = (4000, 2000))
# fig = Figure(resolution = (2000, 3000))
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
# region_df = DataFrames.combine(gb,
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
    pdf = filter(r->r["region"] == major, titration_plot_df)
    pdf = DataFrames.combine(groupby(pdf, [:subregion,:period]),
        Symbol("Δf/f")=>quick_mean_sem=>AsTable, keepkeys=true)
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
    global grid = aog.draw!(fig[i,j], p, axis=(xticklabelrotation = pi/6,))
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

##
tdf = DataFrames.combine(groupby(titration_plot_df, [:subregion,:period]),
           Symbol("Δf/f")=>quick_mean_sem=>AsTable, keepkeys=true)

####################################
## Titration select regions
####################################
W = 183u"mm" / 1.5
# W = 183u"mm" / 2
H = 57u"mm"
dpi = 300 # makie default; we scale up on save
mm2px = mm->Int(round(uconvert(u"inch", mm) * dpi / 1u"inch"))
@show W = mm2px(W)
@show H = mm2px(H)
fig = Figure(resolution=(W, H))

pdf = DataFrames.combine(groupby(titration_plot_df, [:region,:period, :subregion]),
    Symbol("Δf/f")=>quick_mean_sem=>AsTable)

pdf[pdf.subregion .== "Tectum Stratum Periventriculare",:subregion] .= "Optic Tectum"
pdf[pdf.subregion .== "Glyt2 Stripe 3",:subregion] .= "Hindbrain Glyt2 Stripe 3"
pdf[pdf.subregion .== "Gad1b Cluster 16",:subregion] .= "Hindbrain Gad1b Cluster 16"
pdf[pdf.subregion .== "Gad1b Cluster 2",:subregion] .= "Hindbrain Gad1b Cluster 2"
pdf.subregion = CategoricalArray(pdf.subregion, ordered=true)
pdf = filter(row -> row.subregion in REGION_LIST, pdf)
levels!(pdf.subregion, REGION_LIST)
@warn "we initially chose subregions so no collisions; could erroneously combine subregions if violated"
# pdf.period .= replace.(Array(pdf.period), "early" => "1-6")
# pdf.period .= replace.(Array(pdf.period), "mid" => "7-11")
# pdf.period .= replace.(Array(pdf.period), "late" => "12-16")
rename!(pdf, "period" => "# bilateral\ncell stim")
period = "# bilateral\ncell stim"
p2 = Data(pdf) * visual(CrossBar, width=0.5) *
    mapping(:subregion, "Δf/f", "minus_sem", "plus_sem", color=period, dodge=period)

colors = colormap("Blues",7)[2:5]
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

###############################
## save dataframes for plotting in R
###############################

# lol not yet...



###############################
## titration F-test
###############################
periods = ["1", "2-6", "7-11", "12-16"]
region_names = unique(titration_plot_df[!,:fullname])
pvals = map(region_names) do region
    groups = [Float64.(filter(r->(string(r[:period]) == p) && (r.fullname == region), titration_plot_df)[:,"Δf/f"])
        for p in periods]
    # KW = KruskalWallisTest(groups...)
    # pvalue(KW)
    @assert length(groups) == 4
    df1 = DataFrame(dose = 1, score=groups[1])
    df2 = DataFrame(dose = 4, score=groups[2])
    df3 = DataFrame(dose = 9, score=groups[3])
    df4 = DataFrame(dose = 14, score=groups[4])
    # df1 = DataFrame(dose = 1, score=groups[1])
    # df2 = DataFrame(dose = 2, score=groups[2])
    # df3 = DataFrame(dose = 3, score=groups[3])
    # df4 = DataFrame(dose = 4, score=groups[4])
    # not tested; make categorical..?
    # df3[!,:dose] = string.(df3[!,:dose])
    glm_df = vcat(df1, df2, df3, df4)


    nullModel = lm(@formula(score ~ 1), glm_df);
    doseModel = lm(@formula(score ~ 1 + dose), glm_df);

    F = ftest(nullModel.model, doseModel.model);
    F.pval[2]
end
##
# corrected_pvals = MultipleTesting.adjust(pvals, Bonferroni())
corrected_pvals = MultipleTesting.adjust(pvals, BenjaminiHochbergAdaptive())

significant5 = corrected_pvals .< 0.05
significant1 = corrected_pvals .< 0.01
# bonferoni_p5 = 0.05 / length(pvals)
# bonferoni_p1 = 0.01 / length(pvals)
# significant5 = pvals .< bonferoni_p5
# significant1 = pvals .< bonferoni_p1

##
println("=====Significant at p=0.05: $(sum(significant5))/$(length(pvals))=====")
for n in sort(setdiff(region_names[significant5], region_names[significant1])); println(n); end
println("=====Significant at p=0.01: $(sum(significant1))/$(length(pvals))=====")
for n in sort(region_names[significant1]); println(n); end

# raphe p-vals
raphe_idxs = map(x->occursin("Raphe",x), region_names)
collect(zip(region_names[raphe_idxs], corrected_pvals[raphe_idxs]))


# see output at bottom

###############################
## second-order paired Mann–Whitney U test
###############################
second_order_region_names = unique(second_order_plot_df[!,:fullname])
second_order_pvals = map(second_order_region_names) do region
    groups = [Float64.(filter(r->(string(r["stim target"]) == s) &&
        (r.fullname == region), second_order_plot_df)[:,"Δf/f"])
        for s in levels(second_order_plot_df[!,"stim target"])]
    if region == "Rhombencephalon - Gad1b Cluster 2"
        @show groups
    end
    # KW = KruskalWallisTest(groups...)
    # pvalue(KW)
    MWs = [MannWhitneyUTest(groups[s1],groups[s2]) for (s1,s2) in Iterators.product(1:length(groups),1:length(groups))]
    pvals = map(pvalue, MWs)
end
idx = findfirst(second_order_region_names .== "Rhombencephalon - Gad1b Cluster 2")
correct_mat = second_order_pvals[idx]
second_order_pvals_reshape = permutedims(cat(second_order_pvals...,dims=3),(3,1,2))
second_order_pvals_reshape[idx,:,:]
thedims = size(second_order_pvals_reshape)

##
# corrected_pvals = MultipleTesting.adjust(pvals, Bonferroni())
corrected_second_order_pvals = MultipleTesting.adjust(copy(second_order_pvals_reshape[:]), BenjaminiHochbergAdaptive())
corrected_second_order_pvals = reshape(corrected_second_order_pvals, thedims)
corrected_second_order_pvals[idx,:,:]
levels(stim_targets)
##
second_order_significant5 = corrected_second_order_pvals .< 0.05
second_order_significant1 = corrected_second_order_pvals .< 0.01
# bonferoni_p5 = 0.05 / length(pvals)
# bonferoni_p1 = 0.01 / length(pvals)
# significant5 = pvals .< bonferoni_p5
# significant1 = pvals .< bonferoni_p1


# println("=====Significant at p=0.05: $(sum(second_order_significant5))/$(length(corrected_second_order_pvals))=====")
# for n in sort(setdiff(second_order_region_names[second_order_significant5], second_order_region_names[second_order_significant1])); println(n); end
# println("=====Significant at p=0.01: $(sum(second_order_significant1))/$(length(corrected_second_order_pvals))=====")
# for n in sort(second_order_region_names[second_order_significant1]); println(n); end

##

for region in REGION_LIST
    println("======$region======")
    reg = replace(region, "Hindbrain "=>"")
    # idx = findfirst(map(r->occursin(reg,r), second_order_region_names))
    
    # idx = findfirst(map(r -> isnothing(match(r".* Gad1b Cluster 2$", r)),
    #     second_order_region_names))
    # idx = findfirst(map(r -> ~isnothing(match(".* $(region)" * r"$" , r)),
    idx = findfirst(map(r -> ~isnothing(match(Regex(".* $(region)\$") , r)),
        second_order_region_names))

    @show idx
    println(corrected_second_order_pvals[idx,:,:])
end

## for 2021-06-08 titration fish
=====Significant at p=0.05: 60/159=====
Diencephalon -
Diencephalon - Anterior pretectum cluster of vmat2 Neurons
Diencephalon - Dorsal Habenula
Diencephalon - Dorsal Thalamus
Diencephalon - Habenula
Diencephalon - Olig2 Band
Diencephalon - Olig2 Band 2
Diencephalon - Pineal
Diencephalon - Pineal Vmat2 cluster 
Diencephalon - Pretectal Gad1b Cluster
Diencephalon - Pretectal dopaminergic cluster
Diencephalon - Pretectum
Diencephalon - Ventral Habenula
Diencephalon - Ventral Thalamus
Ganglia - Eyes
Ganglia - Lateral Line Neuromast SO3
Mesencephalon -
Mesencephalon - Medial Tectal Band
Mesencephalon - Oculomotor Nucleus nIII
Mesencephalon - Otpb Cluster
Mesencephalon - Oxtl Cluster Sparse
Mesencephalon - Ptf1a Cluster
Mesencephalon - Retinal Arborization Field 7 (AF7)
Mesencephalon - Retinal Arborization Field 9 (AF9)
Mesencephalon - Sparse 6.7FRhcrtR cluster
Mesencephalon - Tecum Neuropil
Mesencephalon - Tegmentum
Mesencephalon - Torus Longitudinalis
Mesencephalon - Vglut2 cluster 1
Rhombencephalon - Anterior Cluster of nV Trigeminal Motorneurons
Rhombencephalon - Caudal Ventral Cluster Labelled by Spinal Backfills
Rhombencephalon - Cerebellum
Rhombencephalon - Cerebellum Gad1b Enriched Areas
Rhombencephalon - Corpus Cerebelli
Rhombencephalon - Gad1b Cluster 1
Rhombencephalon - Gad1b Cluster 16
Rhombencephalon - Gad1b Cluster 17
Rhombencephalon - Gad1b Cluster 2
Rhombencephalon - Gad1b Cluster 5
Rhombencephalon - Glyt2 Cluster 5
Rhombencephalon - Isl1 Cluster 1
Rhombencephalon - Isl1 Cluster 2
Rhombencephalon - Lobus caudalis cerebelli
Rhombencephalon - Neuropil Region 6
Rhombencephalon - Olig2 enriched areas in cerebellum
Rhombencephalon - Otpb Cluster 4
Rhombencephalon - Raphe - Inferior
Rhombencephalon - Raphe - Superior
Rhombencephalon - Rhombomere 1
Rhombencephalon - VII' Facial Motor and octavolateralis efferent neurons
Rhombencephalon - Vglut2 cluster 1
Rhombencephalon - Vglut2 cluster 3
Rhombencephalon - Vmat2 Cluster 3
Telencephalon -
Telencephalon - Pallium
Telencephalon - Subpallial Gad1b cluster
Telencephalon - Subpallial dopaminergic cluster
Telencephalon - Subpallium
Telencephalon - Vglut2 rind
Telencephalon - Vmat2 cluster
=====Significant at p=0.01: 45/159=====
Diencephalon -
Diencephalon - Anterior pretectum cluster of vmat2 Neurons
Diencephalon - Dorsal Habenula
Diencephalon - Dorsal Thalamus
Diencephalon - Habenula
Diencephalon - Olig2 Band
Diencephalon - Olig2 Band 2
Diencephalon - Pineal
Diencephalon - Pineal Vmat2 cluster 
Diencephalon - Pretectal Gad1b Cluster
Diencephalon - Pretectal dopaminergic cluster
Diencephalon - Pretectum
Diencephalon - Ventral Habenula
Diencephalon - Ventral Thalamus
Mesencephalon -
Mesencephalon - Medial Tectal Band
Mesencephalon - Oxtl Cluster Sparse
Mesencephalon - Ptf1a Cluster
Mesencephalon - Retinal Arborization Field 9 (AF9)
Mesencephalon - Tecum Neuropil
Mesencephalon - Tegmentum
Mesencephalon - Torus Longitudinalis
Mesencephalon - Vglut2 cluster 1
Rhombencephalon - Anterior Cluster of nV Trigeminal Motorneurons
Rhombencephalon - Caudal Ventral Cluster Labelled by Spinal Backfills
Rhombencephalon - Cerebellum Gad1b Enriched Areas
Rhombencephalon - Gad1b Cluster 1
Rhombencephalon - Gad1b Cluster 17
Rhombencephalon - Gad1b Cluster 5
Rhombencephalon - Isl1 Cluster 1
Rhombencephalon - Isl1 Cluster 2
Rhombencephalon - Lobus caudalis cerebelli
Rhombencephalon - Neuropil Region 6
Rhombencephalon - Olig2 enriched areas in cerebellum
Rhombencephalon - Otpb Cluster 4
Rhombencephalon - Raphe - Superior
Rhombencephalon - VII' Facial Motor and octavolateralis efferent neurons
Rhombencephalon - Vglut2 cluster 1
Rhombencephalon - Vglut2 cluster 3
Telencephalon -
Telencephalon - Pallium
Telencephalon - Subpallial Gad1b cluster
Telencephalon - Subpallium
Telencephalon - Vglut2 rind
Telencephalon - Vmat2 cluster