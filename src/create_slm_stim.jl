# %% 
# y+1
# %%
using MAT, Dates, DataStructures

# %% DRIVERS

out_name = "stim8-02262020-1249"
# out_name = "registration"
# out_name = "registration_test.txt"

# not positive about z units...
z_offset = 0.
# ntargets = 4
# targets = zeros(ntargets,3) # x (px) by y (px) by z (10^-6 meters aka microns)
# targets = [128. 128. z_offset; 384. 384. z_offset; 128. 384. z_offset; 384. 128. z_offset; 256. 256. 0.]
# targets = [128. 128. z_offset; 384. 384. z_offset; 128. 384. z_offset; 384. 128. z_offset]
targets = [ 369. 272. z_offset
          ; 299. 231. z_offset
          ; 301. 331. z_offset
          ; 360. 274. z_offset
          ; 367. 296. z_offset
          ; 314. 321. z_offset
          ; 250. 366. z_offset
          ; 272. 217. z_offset
          ]
ntargets = size(targets,1)
# %%
# examples
# file = matopen("Z:\\mSLM_B115\\SetupFiles\\CalibrationMasks\\NewMaskFormat\\SLM1_registrationTargets_001_sequential5.mat")
# mat1 = read(file, "experimentS")
# file = matopen("Z:\\mSLM_B115\\SetupFiles\\CalibrationMasks\\NewMaskFormat\\SLM1_registrationTargets_001_sequential5.mat")
# mat2 = read(file, "experimentS")

# %% 
outdir = "O:\\users\\tyler\\slm\\masks"
today_str = Dates.format(Date(now()), DateFormat("Y-mm-dd"))
today_dir = mkpath(joinpath(outdir, today_str))
outfile = joinpath(today_dir, "$out_name.mat")

# %%
stim_struct = OrderedDict{String, Any}()
# field order: im, targets, zrange, targetsGroup
exp_key = "experimentS"
# exp_key = "etb"
stim_struct[exp_key] = Dict{String, Any}()
out_mat = stim_struct[exp_key]

out_mat["im"] = Array{Any}(nothing,1,2)
out_mat["im"][1] = zeros(512,512)
out_mat["im"][2] = zeros(512,512)

out_mat["targets"] = Array{Any}(nothing,1,1)
out_mat["targets"][1] = targets
out_mat["zrange"] = [minimum(targets[:,3]) maximum(targets[:,3])]
# channel 1 (red) and 2 (green), used for GUI preview


# targetsGroup must alternate between Float32(1,3) and zeros(0,0), and end with 2 * zeros(0,0)
out_mat["targetsGroup"] = Array{Any}(nothing,1,ntargets*2+2)
for (i, v) in enumerate(eachrow(targets))
    # SubArrays cannot be written, so we convert to Array
    out_mat["targetsGroup"][i*2-1] = convert(Array,v')
    out_mat["targetsGroup"][i*2] = zeros(0,0)
end
out_mat["targetsGroup"][end-1] = zeros(0,0)
out_mat["targetsGroup"][end] = zeros(0,0)

matwrite(outfile, stim_struct)

out_path = joinpath(today_dir, "$out_name.txt")
txt_io = open(out_path, "w")
print(txt_io, outfile*"\t1")
println("Created $out_path")
close(txt_io)