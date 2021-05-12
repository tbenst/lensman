# lensman
Wielding the power of the lens in galaxies full of neurons

## Getting started
- Install julia
```bash
# (assuming ssh key added to github account)
$ git clone git@github.com:tbenst/lensman.git`
$ cd lensman
$ ml load gurobi # for sherlock, optional
$ julia
julia>ENV["GUROBI_HOME"] = "/opt/gurobi911" # optional
julia>ENV["GUROBI_HOME"] = "/share/software/user/restricted/gurobi/9.0.3" # sherlock
julia>ENV["PYTHON"] = "/home/tyler/opt/anaconda3/envs/jl-1.6/bin" # optional
julia>ENV["PYTHON"] = "/home/groups/deissero/tyler/anaconda/bin/python" # optional
julia>]
pkg>registry add https://github.com/HolyLab/HolyLabRegistry.git
pkg>activate .
pkg>instantiate
```

For Julia 1.6, need to use this workaround for MKL: https://github.com/JuliaLinearAlgebra/MKL.jl/issues/60#issuecomment-808329869
### auto activate environment
for convenience, you can skip `activate .` by creating `$HOME/.julia/config/startup.jl`:
```julia
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
```

## adding registries / unregistered packages (if needed)
```julia
julia>]
pkg> activate .
pkg> registry add https://github.com/HolyLab/HolyLabRegistry.git
pkg> add https://github.com/timholy/ANTsRegistration.jl
pkg> add https://github.com/JuliaFolds/GeneratorsX.jl
```

## fix old glibc issue
https://github.com/JuliaGL/GLFW.jl/issues/198#issue-574359095