ENV["DISPLAY"] = "localhost:12"
##
using AlgebraOfGraphics, CairoMakie, PalmerPenguins, DataFrames

set_aog_theme!() #src

penguins = dropmissing(DataFrame(PalmerPenguins.load()))

p1 = data(penguins) * visual(Violin) *
    mapping(:species, :bill_depth_mm, color=:sex, dodge=:sex)
draw(p1, axis=(xscale=log,))
