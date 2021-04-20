using Lensman
using Unitful: μm, m, s, uconvert

## SLM 1
offset2Microns = 120μm
slmNum = 1
offset = float(uconvert(m, offset2Microns)) / m
targets1 = [128. 128. offset; 384. 384. offset; 128. 384. offset; 384. 128. offset;]

offsetStr = replace(string(offset2Microns), " μm"=>"_um")
create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2021-02-01_4square_SLM$(slmNum)_offset_$(offsetStr)",
    slmNum=slmNum)


## SLM 2
offset2Microns = 0μm
slmNum = 2
offset = float(uconvert(m, offset2Microns)) / m
x1,x2,y1,y2 = 128., 384., 128., 384.
# x1,x2,y1,y2 = 10., 502., 10., 502.
# x1,x2,y1,y2 = 50., 512. -50, 50., 512 - 50.
targets = [[x1 y1 offset], [x2 y2 offset], [x1 y2 offset], [x2 y1 offset]]
xy_offset = 50
# targets_center_list = [[x1+xy_offset y1+xy_offset], [x2-xy_offset y2-xy_offset], [x1+xy_offset y2-xy_offset], [x2-xy_offset y1+xy_offset]]
targets_center_list = nothing

offsetStr = replace(string(offset2Microns), " μm"=>"_um")
create_slm_stim(targets,
    "/mnt/deissero/users/tyler/slm/masks/2021-03-30_4square_1concurrent_SLM$(slmNum)",
    slmNum=slmNum, targets_center_list=targets_center_list)

## SLM 2 (galvo pan)
offset2Microns = 0μm
slmNum = 2
offset = float(uconvert(m, offset2Microns)) / m
# x1,x2,y1,y2 = 128., 384., 128., 384.
# x1,x2,y1,y2 = 10., 502., 10., 502.
x1,x2,y1,y2 = 50., 512. -50, 50., 512 - 50.
targets = [[x1 y1 offset], [x2 y2 offset], [x1 y2 offset], [x2 y1 offset]]
xy_offset = 30
# targets_center_list = [[x1+xy_offset y1+xy_offset], [x2-xy_offset y2-xy_offset], [x1+xy_offset y2-xy_offset], [x2-xy_offset y1+xy_offset]]
# only plus
targets_center_list = [[x1+xy_offset y1+xy_offset], [x2+xy_offset y2+xy_offset], [x1+xy_offset y2+xy_offset], [x2-xy_offset y1+xy_offset]]

offsetStr = replace(string(offset2Microns), " μm"=>"_um")
create_slm_stim(targets,
    "/mnt/deissero/users/tyler/slm/masks/2021-03-30_4corner_galvo_pan_SLM$(slmNum)_xyoffset_+$(xy_offset)",
    slmNum=slmNum, targets_center_list=targets_center_list)


## SLM 2 (galvo pan) with pattern (1cell stim)
offset2Microns = 0μm
slmNum = 2
offset = float(uconvert(m, offset2Microns)) / m
# x1,x2,y1,y2 = 128., 384., 128., 384.
# x1,x2,y1,y2 = 10., 502., 10., 502.
x1,x2,y1,y2 = 50., 512. -50, 50., 512 - 50.
targets = [[x1 y1 offset], [x2 y2 offset], [x1 y2 offset], [x2 y1 offset]]
xy_offset = 30
# targets_center_list = [[x1+xy_offset y1+xy_offset], [x2-xy_offset y2-xy_offset], [x1+xy_offset y2-xy_offset], [x2-xy_offset y1+xy_offset]]
# only plus
targets_center_list = [[x1+xy_offset y1+xy_offset], [x2+xy_offset y2+xy_offset], [x1+xy_offset y2+xy_offset], [x2-xy_offset y1+xy_offset]]

offsetStr = replace(string(offset2Microns), " μm"=>"_um")
create_slm_stim(targets,
    "/mnt/deissero/users/tyler/slm/masks/2021-03-30_4corner_galvo_pan_SLM$(slmNum)_xyoffset_+$(xy_offset)",
    slmNum=slmNum, targets_center_list=targets_center_list)



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
slmNum = 2
topOffset = float(uconvert(m, -50μm)) / m
bottomOffset = float(uconvert(m, 50μm)) / m
targets1 = [128. 128. topOffset; 384. 384. topOffset; 128. 384. topOffset; 384. 128. topOffset;]
targets2 = [128. 256. bottomOffset; 256. 128. bottomOffset; 256. 384. bottomOffset; 384. 256. bottomOffset;]
topOffsetStr = replace(string(topOffset), " μm"=>"_um")
bottomOffsetStr = replace(string(bottomOffset), " μm"=>"_um")
create_slm_stim([targets1, targets2],
    "/mnt/deissero/users/tyler/slm/masks/2021-02-15_8square_SLM$(slmNum)_offset_$(bottomOffsetStr)_$(topOffsetStr)",
    slmNum=slmNum)

##
topOffset = float(uconvert(m, 20μm)) / m
bottomOffset = float(uconvert(m, 60μm)) / m
targets1 = [128. 128. topOffset; 384. 384. topOffset; 128. 384. topOffset; 384. 128. topOffset;]
targets2 = [128. 256. bottomOffset; 256. 128. bottomOffset; 256. 384. bottomOffset; 384. 256. bottomOffset;]

create_slm_stim([targets1, targets2],
    "/mnt/deissero/users/tyler/slm/masks/2020-11-16_registration_8square_20_60um")
##
zOffset = float(uconvert(m, 10μm)) / m
slmNum = 2
first = round(typeof(1μm), uconvert(μm, zOffset*(1-8)m)) / μm
last = round(typeof(1μm), uconvert(μm, zOffset*(16-8)m)) / μm
targets1 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 1:4]...)')
targets2 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 5:8]...)')
targets3 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 9:12]...)')
targets4 = copy(hcat([[256, 16+32*(i-1), zOffset*(i-8)] for i in 13:16]...)')

create_slm_stim([targets1, targets2, targets3, targets4],
    "/mnt/deissero/users/tyler/slm/masks/2021-02-09_SLM$(slmNum)_12line_zsweep_$(first)_$(last)",
    slmNum=slmNum)