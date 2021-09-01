ENV["DISPLAY"] = "localhost:11.0"

##
using DataFrames, Statistics, Arrow, AlgebraOfGraphics, CairoMakie
Data = AlgebraOfGraphics.data

##
arrow_path = "/scratch/allan/b115/2021-06-08_rsChRmine_h2b6s/fish2/h5_output/TSeries-lrhab-titration-1232021-06-21_6pm.ty.h5/PerVoxelLSTM_actually_shared-separate_bias_hidden-2021-06-21_6pm/regions_df.arrow"
regions_df = DataFrame(Arrow.Table(arrow_path))

##
mask_names = [
    "Diencephalon - Habenula",
    "Mesencephalon - Tectum Stratum Periventriculare",  # only one
    "Rhombencephalon - Cerebellum",  # choosing most general one
    # "Rhombencephalon - Raphe - Inferior",
    "Rhombencephalon - Raphe - Superior",
    "Rhombencephalon - Rhombomere 1"
]
plot_df = filter(row -> any(mask_names .== row.region), regions_df)
unique(plot_df[:, :region])

## 
fig = Figure(resolution=(1200, 500))
p1 = Data(filter(row -> row.stim == "1" || row.stim == "16", plot_df)) * visual(BoxPlot, show_notch=true) *
    mapping(:region, "Δf/f", color=:hemisphere, dodge=:hemisphere, col=:stim)
grid = draw!(fig, p1, axis=(xticklabelrotation = pi / 6,))
legend!(fig[1,end + 1], grid)
fig
##
fig = Figure(resolution=(1200, 500))
df = combine(groupby(plot_df, [:region, :stim]), "Δf/f" => mean)
df[:,:stim_group] = parse.(Int, df[:,:stim])
p1 = Data(df) * visual(Scatter) *
    mapping(:stim_group, "Δf/f_mean", color=:region)
grid = draw!(fig, p1, axis=(xticklabelrotation = pi / 6,))
legend!(fig[1,end + 1], grid)
fig


# First do *no* grouping and just define set of 
