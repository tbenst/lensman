# lensman
Wielding the power of the lens in galaxies full of neurons

## Getting started
- Install julia
```bash
# (assuming ssh key added to github account)
$ git clone git@github.com:tbenst/lensman.git`
$ cd lensman
$ julia
julia>]
pkg>activate .
pkg>instantiate
```

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