using Dates
using Unitful: μm, m, s, mW

function getzoffset(exp_date, slmNum)
    if slmNum==1
        if exp_date < Date(2021,1,11)
            offset = float(uconvert(m, 45μm)) / m # prior to 2021 / starting on ...12/15...? should check...
        elseif exp_date < Date(2021,2,2)
            offset = float(uconvert(m, 48μm)) / m # since 2020-01-11
        else 
            offset = float(uconvert(m, 45μm)) / m
        end
    elseif slmNum==2
        offset = float(uconvert(m, 0μm)) / m
    end
    zOffset = offset * 1e6
    return zOffset
end

"TODO: could potentially read this from sean's files...?"
function spiral_size(exp_date, microscope_lateral_unit)
    if exp_date < Date(2020,12,4)
        @warn "Need to double check spiral size for old experiments"
        15μm * (14.4/25) / microscope_lateral_unit
    elseif exp_date < Date(2021,3,9)
        7μm * (14.4/25) / microscope_lateral_unit
    else
        7μm / microscope_lateral_unit
    end
end

# TODO: read objective from Prarie view xml
# TODO: measure size in Prarie view
function microscope_lateral_unit(W)
    @info "Assume 25x objective & 1x zoom"
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