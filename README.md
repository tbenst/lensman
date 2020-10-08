# lensman
Wielding the power of the lens in galaxies full of neurons

## TO RUN
julia --project=. .\create_slm_stim.jl

## Install
- add repository:
- 
```julia
julia>]
pkg> registry add https://github.com/HolyLab/HolyLabRegistry.git
pkg> activate .
```


### nixos
- enter directory to launch direnv shell (lorri)
- `code .`
- that's it!

downside: this won't work with remote-ssh :/

# registration
- https://github.com/timholy/ANTsRegistration.jl


## papers:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5597853/#!po=16.6667



### ANTs debugging....
#### example command
/home/tyler/.julia/packages/ANTsRegistration/GZjzQ/deps/usr/bin/antsRegistration -d 3 --transform 'Affine[0.1]' --metric 'MI[/tmp/tyler/ANTs/3muoxnnd66.nrrd,/tmp/tyler/ANTs/dTmlyz29Lu.nrrd,1.0,32,Regular,0.25]' --convergence '[1000x500x250x0,1.0e-6,10]' --shrink-factors 12x8x4x2 --smoothing-sigmas 4x3x2x1vox --output '[/tmp/tyler/ANTs/ztW88wKbBF_warp,/tmp/tyler/ANTs/ztW88wKbBF.nrrd]'

LOOK for ztW88wKbBF.nii.gz. Then apply with `antsApplyTransforms`

-j, --save-state saveSateAsTransform
          A boolean specifying whether or not the composite transform resulted from the 
          current state of the registration should be written to an hdf5 composite file. 
          It is specially usefull if we want to save the current state of a SyN 
          registration to the disk, so we can load and restore that later to continue the 
          next registration directly started from the same state that it was stored. The 
          output of this flag will be the same as the write-composite-transform, unless 
          the last transform is a SyN transform. In that case, the inverse displacement 
          field of the SyN transform is also added to the output composite transform. To 
          load the written file this flag, we must use restore-state.


## NixOS (IMPORTANT)

Ants must be installed in develop mode, then after bulid, modify deps.jl (/home/tyler/.julia/dev/ANTsRegistration/deps/deps.jl) 
```
const ants_dir = read(`nix eval --raw "nixpkgs.ants"`, String)
const ants = joinpath(ants_dir, "bin/ANTS")
const antsRegistration = joinpath(ants_dir, "bin/antsRegistration")
const antsMotionCorr = joinpath(ants_dir, "bin/antsMotionCorr")
```