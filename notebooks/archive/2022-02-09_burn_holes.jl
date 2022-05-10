using Lensman
using Unitful: μm, m, s, uconvert

# can't burn 30um down or -30um up

## SLM 1 (works in B118)
offset2Microns = 10μm
slmNum = 1
offset = float(uconvert(m, offset2Microns)) / m
targets1 = [128.0 128.0 offset; 384.0 384.0 offset; 128.0 384.0 offset; 384.0 128.0 offset]

offsetStr = replace(string(offset2Microns), " μm" => "_um")
create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2022-02-09_4square_SLM$(slmNum)_offset_$(offsetStr)",
    slmNum=slmNum)

##
zOffset = float(uconvert(m, 10μm)) / m
slmNum = 1
zfirst = round(typeof(1μm), uconvert(μm, zOffset * (1 - 8)m)) / μm
zlast = round(typeof(1μm), uconvert(μm, zOffset * (16 - 8)m)) / μm
x = 128
targets1 = copy(hcat([[x, 16 + 32 * (i - 1), zOffset * (i - 8)] for i in 1:4]...)')
targets2 = copy(hcat([[x, 16 + 32 * (i - 1), zOffset * (i - 8)] for i in 5:8]...)')
targets3 = copy(hcat([[x, 16 + 32 * (i - 1), zOffset * (i - 8)] for i in 9:12]...)')
targets4 = copy(hcat([[x, 16 + 32 * (i - 1), zOffset * (i - 8)] for i in 13:16]...)')

create_slm_stim([targets1, targets2, targets3, targets4],
    # "/mnt/deissero/users/tyler/slm/masks/2022-02-09_SLM$(slmNum)_12line_zsweep_$(zfirst)_$(zlast)_x=$(x)",
    "/mnt/deissero/users/tyler/slm/masks/debug4",
    slmNum=slmNum)

create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/debug1",
    slmNum=slmNum)
create_slm_stim([targets2],
    "/mnt/deissero/users/tyler/slm/masks/debug2",
    slmNum=slmNum)
create_slm_stim([targets3],
    "/mnt/deissero/users/tyler/slm/masks/debug3",
    slmNum=slmNum)
create_slm_stim([targets4],
    "/mnt/deissero/users/tyler/slm/masks/debug4",
    slmNum=slmNum)


##
# fails!
t1 = [128.0 16.0 0.0
    128.0 48.0 0.0
    128.0 80.0 -1.0e-5
    128.0 112.0 -2.0e-5
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d1a",
    slmNum=slmNum)
t1 = [128.0 16.0 0.0
    128.0 48.0 0.0
    128.0 80.0 -1.0e-5
    128.0 112.0 -3.0e-5
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d2b",
    slmNum=slmNum)


t1 = [128.0 16.0 0.0
    128.0 48.0 0.0
    128.0 80.0 1.0e-5
    128.0 112.0 2.0e-5
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d2a",
    slmNum=slmNum)
t1 = [128.0 16.0 0.0
    128.0 48.0 0.0
    128.0 80.0 1.0e-5
    128.0 112.0 3.0e-5
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d2b",
    slmNum=slmNum)


##

# fail
t1 = [128.0 16.0 -4.0e-5
    128.0 48.0 -3.0e-5
    128.0 80.0 -2.0e-5
    128.0 112.0 -1.0e-5
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d3",
    slmNum=slmNum)

# okay 
t1 = [128.0 16.0 0.0
    128.0 48.0 0.0
    128.0 80.0 -1.0e-5
    128.0 112.0 -1.0e-5
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d4",
    slmNum=slmNum)

t1 = [128.0 16.0 0.0
    128.0 48.0 0.0
    128.0 80.0 0.0
    128.0 112.0 0.0
]

create_slm_stim([t1],
    "/mnt/deissero/users/tyler/slm/masks/d2",
    slmNum=slmNum)



##
create_slm_stim([targets3],
    "/mnt/deissero/users/tyler/slm/masks/2022-02-09_SLM$(slmNum)_4line_zsweep_$(zfirst)_$(zlast)_x=$(x)",
    slmNum=slmNum)

## succeeds
offset = 0.0
targets1 = [128.0 272.0 offset; 304.0 384.0 offset; 128.0 336.0 offset; 384.0 368.0 offset]
targets2 = [136.0 260.0 offset; 304.0 350.0 offset; 136.0 350.0 offset; 384.0 350.0 offset]
create_slm_stim([targets1, targets2],
    "/mnt/deissero/users/tyler/slm/masks/2022-02-09_debug3",
    slmNum=slmNum)

create_slm_stim([targets1],
    "/mnt/deissero/users/tyler/slm/masks/2022-02-09_debug1",
    slmNum=slmNum)

create_slm_stim([targets2],
    "/mnt/deissero/users/tyler/slm/masks/2022-02-09_debug2",
    slmNum=slmNum)

## works
offset = 0.0
targets1 = [110.0 110.0 offset; 120.0 120.0 offset; 130.0 130.0 offset; 140.0 140.0 offset]
targets2 = [210.0 210.0 offset; 220.0 220.0 offset; 230.0 230.0 offset; 240.0 240.0 offset]
targets3 = [310.0 310.0 offset; 320.0 320.0 offset; 330.0 330.0 offset; 340.0 340.0 offset]
targets4 = [410.0 410.0 offset; 420.0 420.0 offset; 430.0 430.0 offset; 440.0 440.0 offset]
create_slm_stim([targets1, targets2, targets3, targets4],
    "/mnt/deissero/users/tyler/slm/masks/2022-02-09_debug3",
    slmNum=slmNum)



##
offset = 0.0
slmNum = 1
targets = []

grid_space = 512 / 9
xs = Float64.(collect(32:64:512))
ys = Float64.(collect(32:64:512))
for y in ys
    push!(targets, copy(hcat([[x, y, 0.0] for x in xs[1:4]]...)'))
    push!(targets, copy(hcat([[x, y, 0.0] for x in xs[5:8]]...)'))
end

create_slm_stim(targets,
    "/mnt/deissero/users/tyler/slm/masks/2022-02-22_8x8_grid",
    slmNum=slmNum)

##
offset = 0.0
slmNum = 1
targets = []

grid_space = 512 / 9
xs = Float64.(collect(32:64:512)) .+ rand(collect(-3:3), 8)
ys = Float64.(collect(32:64:512)) .+ rand(collect(-3:3), 8)
for y in ys
    push!(targets, copy(hcat([[x, y, 0.0] for x in xs[1:4]]...)'))
    push!(targets, copy(hcat([[x, y, 0.0] for x in xs[5:8]]...)'))
end
targets
create_slm_stim(targets,
    "/mnt/deissero/users/tyler/slm/masks/2022-05-04_8x8_dithered",
    slmNum=slmNum)

##
offset = 0.0
slmNum = 1
targets = []

grid_space = 512 / 9
xs = Float64.(collect(32:64:512)) .+ rand(collect(-3:3), 8)
ys = Float64.(collect(32:64:512)) .+ rand(collect(-3:3), 8)
for x in xs
    push!(targets, copy(hcat([[x, y, 0.0] for y in ys[1:4]]...)'))
    push!(targets, copy(hcat([[x, y, 0.0] for y in ys[5:8]]...)'))
end
targets
create_slm_stim(targets,
    "/mnt/deissero/users/tyler/slm/masks/2022-05-04_8x8_dithered_vertical",
    slmNum=slmNum)
