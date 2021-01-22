
## PROTOCOL:
- Olympus 25x 1.05NA, 1.5x zoom, 1024x1024
- Disable micro manipulator!!
- zstack, ~111 slices, 820nm, 64avg, 2um steps, 1x zoom (maybe more slices for 6dpf, maybe 0.8x for better registration since more hindbrain...)
    - 180 to 210 pockels (~20mW-30mW), Bruker default exponential correction
    - takes ~10min
    - 1.5x zoom fails to register to Zbrain (only tested at 11dpf)
- register 
  - takes ~5min
- choose imaging planes w/ code
- take 128avg image @ 1.5x for these planes
- registration part II: get mask for high-res (hopefully just translate / scale?)
- code selects cells inside desired region masks

## LESSONS
- 25x:1x is hard to segment (althoug havent tried latest pipeline...)
- 1.5x is nice to segment when shallow. Hard when at raphe level for 11 dpf
- 1.5x gives very bad zbrain registration
- 1x gives reasonable however
- maybe tray 32 frames for coarse atlas to speed things up..
- 
## TO RUN
julia --project=. .\create_slm_stim.jl

## Install
- add repository:
- 
```julia
julia>]
pkg> activate .
pkg> registry add https://github.com/HolyLab/HolyLabRegistry.git
pkg> add https://github.com/timholy/ANTsRegistration.jl
pkg> add https://github.com/JuliaFolds/GeneratorsX.jl
```


### nixos
- enter directory to launch direnv shell (lorri)
- `code .`
- that's it!

downside: this won't work with remote-ssh :/

### remote-ssh
- need to open separate ssh tunnel with X11 forwarding
- `echo $DISPLAY`
- in julia, `ENV["DISPLAY"] = "localhost:11.0"`


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

### GLFW
GLFW is fixed by specifying "julia.executablePath": "/home/tyler/.nix-profile/bin/julia" !!! (but only locally I think, otherwise need virtualgl)

## Timings
position fish, setup Zseries, phase offset - 8 min
zseries - 9 min
file conversion - 1:30 min
file transfer, get filename into julia, read in - 2 min
ANTs - 8 min
segmentation / choose planes / inspect - 2:30 min
setup z-plane & start Tseries - 2:15 min
TOTAL: ~34 min

biggest optimization: human-selected planes for first Tseries.
Could start first Tseries ~20 min faster, but if Tseries only
takes 10 min then only save 10 min.

other optimizations:
- save 4-6 min with faster CPU / optimize ANTs
- fewer Zseries planes or 32 frame avg = 5 min saving
- start/stream zseries from Julia = 3:30min saving
- automate Zseries construction = 4-5 min saving
  (especially if just put fluorophore at top and go)


### shared drive
sudo mount -t cifs -o username=mslm_cnc115 //171.65.101.227/CNC_B115_mSLM /mnt/b115_mSLM/

## rsync
name="2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/TSeries-1024cell-32concurrent-047"
rsync -avP sher:'/oak/stanford/groups/deissero/users/tyler/b115/2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/*.ty.h5' /data/dlab/b115/2021-01-12_chrmine-kv2.1_h2b6s_7dpf/fish2_chrmine/