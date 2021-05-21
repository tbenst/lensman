import Unitful: MHz, mW, nm, J, cm, Hz, fs, mJ, ns
# conclusion from 2021-05-11: use 2MHz and likely fine until 50 ppc
# rep_rate = 1MHz
# rep_rate = 2MHz
rep_rate = 4MHz
# rep_rate = 10MHz
# avg_power = 340mW
# avg_power = 500mW
avg_power = 300mW
pulse_width = 258fs
power_variance_discount = 0.8

λ = 1035nm
NA = 1.05 * 0.75 # adjust for not filling entire objective...

# https://www.nature.com/articles/s41467-018-08179-6
# 2J/cm^2 at 60 fs in mouse

# we convert damage threshold by pulse width using sqrt approximation
# https://pubmed.ncbi.nlm.nih.gov/2705929/

damage_threshold_upper = 2J/cm^2
threshold_pulse_width = 60fs
adj_damage_threshold_upper = damage_threshold_upper * sqrt(pulse_width)/sqrt(threshold_pulse_width)

# we estimate 1.25J/cm^2 in larval zebrafish
damage_threshold_est = 1.25J/cm^2
threshold_pulse_width = 60fs
adj_damage_threshold_est = damage_threshold_est * sqrt(pulse_width)/sqrt(threshold_pulse_width)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3123732/
# 100fs: 1J/cm^2
# 1ns: 100J/cm^2
# sqrt rule for conversion by pulse width:
@assert uconvert(NoUnits,sqrt(1ns)/sqrt(100fs)) == 100
damage_threshold_lower = 1.0J/cm^2
threshold_pulse_width = 100fs
adj_damage_threshold_lower = damage_threshold_lower * sqrt(pulse_width)/sqrt(threshold_pulse_width)

pulse_energy = avg_power / rep_rate
diffraction_limit = (λ / (2*NA))^2
engery_per_area = pulse_energy / diffraction_limit
max_energy_per_area = uconvert(J/cm^2, engery_per_area)
ppc_lower = power_variance_discount * adj_damage_threshold_lower/max_energy_per_area*1000
ppc_est = power_variance_discount * adj_damage_threshold_est/max_energy_per_area*1000
ppc_upper = power_variance_discount * adj_damage_threshold_upper/max_energy_per_area*1000
println("lower bound safe power per cell: $(ppc_lower)")
println("estimated bound safe power per cell: $(ppc_est)")
println("upper bound safe power per cell: $(ppc_upper)")