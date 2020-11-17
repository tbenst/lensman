using Lensman
using Unitful: μm, m, s, uconvert

## to burn at etl=0 if using calibration circa 2020-11-16, need +30 offset
offset = float(uconvert(m, 30μm)) / m
targets1 = [128. 128. offset; 384. 384. offset; 128. 384. offset; 384. 128. offset;]

create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-16_registration_4square_30um")


##

topOffset = float(uconvert(m, -20μm)) / m
bottomOffset = float(uconvert(m, 40μm)) / m
targets1 = [128. 128. topOffset; 384. 384. topOffset; 128. 384. topOffset; 384. 128. topOffset;]
targets2 = [128. 256. bottomOffset; 256. 128. bottomOffset; 256. 384. bottomOffset; 384. 256. bottomOffset;]

create_slm_stim([targets1, targets2],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-16_registration_8square_-20_40um")

##
topOffset = float(uconvert(m, 20μm)) / m
bottomOffset = float(uconvert(m, 60μm)) / m
targets1 = [128. 128. topOffset; 384. 384. topOffset; 128. 384. topOffset; 384. 128. topOffset;]
targets2 = [128. 256. bottomOffset; 256. 128. bottomOffset; 256. 384. bottomOffset; 384. 256. bottomOffset;]

create_slm_stim([targets1, targets2],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-16_registration_8square_20_60um")
##
zOffset = float(uconvert(m, 10μm)) / m
targets1 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-5)] for i in 1:4]...)')
targets2 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-5)] for i in 5:8]...)')
targets3 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-5)] for i in 9:12]...)')
targets4 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-5)] for i in 13:16]...)')

create_slm_stim([targets1, targets2, targets3, targets4],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-16_registration_12line_zsweep_-40_110")