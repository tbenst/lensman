using Lensman
using Unitful: μm, m, s, uconvert

## SLM 2
offset2Microns = 120μm
slmNum = 1
offset = float(uconvert(m, offset2Microns)) / m
targets1 = [128. 128. offset; 384. 384. offset; 128. 384. offset; 384. 128. offset;]

offsetStr = replace(string(offset2Microns), " μm"=>"_um")
create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2021-02-01_4square_SLM$(slmNum)_offset_$(offsetStr)",
    slmNum=slmNum)


## to burn at etl=0 if using calibration circa 2020-12-15, need +50 offset
# changed to 48mu on 2020-01-11
slmNum = 1
offsetMicrons = 48μm
offset = float(uconvert(m, offsetMicrons)) / m
targets1 = [128. 128. offset; 384. 384. offset; 128. 384. offset; 384. 128. offset;]

offsetStr = replace(string(offsetMicrons), " μm"=>"_um")
create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2021-01-11_4square_SLM$(slmNum)_offset_$(offsetStr)",
    slmNum=slmNum)


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
targets1 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 1:4]...)')
targets2 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 5:8]...)')
targets3 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 9:12]...)')
targets4 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 13:16]...)')

create_slm_stim([targets1, targets2, targets3, targets4],
    "/mnt/deissero/users/tyler/slm/masks/2020-12-15_registration_12line_zsweep_-70_80")