# this notbeook saves a .nrrd file for a given recording:
# e.g. /scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_fixed.nrrd
# for zseries
# and /scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_moving.nrrd
# for h2b_zbrain
using Lensman
import Lensman: @pun
##
resources = Resources()
##
recording1 = Recordings[
        "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-titration-192trial-062"
        # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-titration-192trial-061"
    ](;resources...,
        registration_type = :dont_run
    );
##
@pun (fish_dir, zseries) = recording1;
fixedname = joinpath(fish_dir, "zseries.nrrd")
save(fixedname, zseries)
fixedname