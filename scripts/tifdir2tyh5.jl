import Pkg
Pkg.activate(raw"O:\users\tyler\code\lensman")
using Lensman
# tif_dir = raw"F:\tyler\2021-01-12_chrmine-kv2.1_h2b6s_7dpf\fish2_chrmine\TSeries-1024cell-32concurrent-5ppc-048"

if ~size(ARGS,1)==1
    println(raw"Usage: julia O:\users\tyler\code\lensman\notebooks\2021-01-12_save_tyh5.jl tifdir")
else
    tifDir = ARGS[1]
    tseries = loadTseries(tifDir);
    write_experiment_to_tyh5(tseries, tifDir*".ty.h5", compression_level=0)
end