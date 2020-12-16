@assert size(ARGS,1)==1
using Lensman
# example
# tifDir = "/mnt/deissero/users/tyler/b115/2020-10-28_elavl3-chrmine-Kv2.1_h2b6s_8dpf/fish1/TSeries_lrhab_raphe_40trial_part_2-040"
tifDir = ARGS[1]

slmDir = "/mnt/b115_mSLM/mSLM_B115/SetupFiles/Experiment/" # TODO: automatic backup of dir to OAK...


expName = splitpath(tifDir)[end]
fishDir = joinpath(splitpath(tifDir)[1:end-1]...)
dataFolders = splitpath(tifDir)
xmlPath = joinpath(dataFolders..., dataFolders[end] * ".xml")
expDate, frameRate, etlVals = getExpData(xmlPath)
println(etlVals)