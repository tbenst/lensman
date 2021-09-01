##
ENV["DISPLAY"] = "localhost:12.0"
##
using GLFW; window = GLFW.CreateWindow(640, 480, "GLFW.jl")
##

using GLMakie, AbstractPlotting

## simple scatter
x = 1:10
y = 1:10
sizevec = x ./ 10
scene = scatter(x, y, markersize = sizevec)

## animation

scene = Scene()
mytime = Node(0.0)
f(v, t) = sin(v + t)
g(v, t) = cos(v + t)
lines!(scene,lift(t -> f.(range(0, stop = 2pi, length = 50), t), mytime),color = :blue)
lines!(scene,lift(t -> g.(range(0, stop = 2pi, length = 50), t), mytime),color = :orange)
record(scene, "output.mp4", range(0, stop = 4pi, length = 100)) do i
    mytime[] = i
end
