using Dates
using Unitful: μm, m, s, mW

function getzoffset(exp_date, slmNum)
    if slmNum == 1
        if exp_date < Date(2021, 1, 11)
            offset = float(uconvert(m, 45μm)) / m # prior to 2021 / starting on ...12/15...? should check...
        elseif exp_date < Date(2021, 2, 2)
            offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
        elseif exp_date < Date(2022, 1, 1) # B115
        else
            # B118
            offset = float(uconvert(m, 0μm)) / m
        end
    elseif slmNum == 2
        offset = float(uconvert(m, 0μm)) / m
    end
    zOffset = offset * 1e6
    return zOffset
end

"TODO: could potentially read this from sean's files...?"
function spiral_size(exp_date, microscope_lateral_unit)
    if exp_date < Date(2020, 12, 4)
        @warn "Need to double check spiral size for old experiments"
        15μm * (14.4 / 25) / microscope_lateral_unit
    elseif exp_date < Date(2021, 3, 9)
        7μm * (14.4 / 25) / microscope_lateral_unit
    else
        7μm / microscope_lateral_unit
    end
end

# TODO: read objective from Prarie view xml
# TODO: measure size in Prarie view
function microscope_lateral_unit(W)
    @info "Assume 25x objective & 1x zoom"
    @warn "should use read_microns_per_pixel instead!"
    if W == 1024
        0.6299544139175637μm
    elseif W == 512
        0.6299544139175637μm * 2
    else
        0.6299544139175637μm * 1024 / W
    end
end

function slmpower(exp_date)
    if exp_date < Date("2021-02-24")
        slm1Power = 570mW
        slm2Power = 380mW
    elseif exp_date <= Date("2021-03-01")
        slm1Power = 450mW
        slm2Power = 375mW
    else
        # slm1Power = 377mW
        # slm2Power = 305mW
        slm1Power = 450mW
        slm2Power = 450mW
    end
    slm1Power, slm2Power
end

"""Return min and max voltage for X, Y galvos.

These numbers are used in markpoints xml files. Can identify by taking a
max resolution image, e.g. 1024x1024, and looking at i.e.

```
    <PVStateValue key="maxVoltage">
      <IndexedValue index="XAxis" value="4.33123017" />
      <IndexedValue index="YAxis" value="4.554884304" />
    </PVStateValue>
```
"""
function markpoints_magic_numbers(room::String)
    @assert in(room, ["B113", "B115", "B118"])
if room == "B115"
    minX = -4.17
    minY = -4.62
    maxX = 4.17
    maxY = 4.62
elseif room == "B118"
    # minX = -4.01217945
    # minY = -4.720407024
    # maxX = 4.33123017
    # maxY = 4.554884304
    # post 2022-05-04
    minX = -4.17170481
    minY = -4.17170481
    maxX = 4.17170481
    maxY = 4.17170481
elseif room == "B113"
    minX = -7.6
    maxX = 7.6
    minY = -8.3
    maxY = 8.3
else
    error("bad room string")
end
    (minX, maxX), (minY, maxY)
end