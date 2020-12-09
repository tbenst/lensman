using Compose, Gadfly, Images, TestImages
 
struct BackgroundImage <: Gadfly.GuideElement
    ctx::Context
end

function BackgroundImage(img::AbstractArray{<: Colorant, 2})
    yrange, xrange = Base.axes(img)
    BackgroundImage(compose(context(),
              bitmap("image/png",
                     repr("image/png",img), first(xrange),first(yrange),last(xrange),last(yrange))))
end

function Gadfly.render(guide::BackgroundImage, theme::Gadfly.Theme,
                aes::Gadfly.Aesthetics)
    ctx = compose(context(), svgclass("geometry"), guide.ctx)
    return [Guide.PositionedGuide([ctx], 0, Guide.under_guide_position)]
end

function cartIdx2Array(x::Vector{<:CartesianIndex})
    # make tall matrix
    permutedims(reduce(hcat, collect.(Tuple.(x))), (2,1))
end

img = testimage("lighthouse")
img1 = Gray.(img)
points = cartIdx2Array(findall(img1 .> 0.95))
coords = Coord.cartesian(xmin=0, xmax=size(img1,2), ymin=-size(img,1),ymax=0)

p = plot(x=points[:,2], y=-points[:,1], Geom.point,BackgroundImage(img), coords)
draw(SVG("test.svg"), p)