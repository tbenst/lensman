# this notbeook saves a .nrrd file for a given recording:
# e.g. /scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_fixed.nrrd
# for zseries
# and /scratch/b115/2021-06-08_rsChRmine_h2b6s/fish2/20210826T163344881_moving.nrrd
# for h2b_zbrain
using Lensman, NRRD, FileIO
import Lensman: @pun
##
resources = Resources()
##
recording1 = Recordings[
        # "2021-06-02_rsChRmine-h2b6s/fish2/TSeries-titration-192trial-062"
        # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-titration-192trial-061"
        # "2021-06-01_rsChRmine_h2b6s/fish3/TSeries-titration-192trial-061"
        # "2021-06-15_rsChRmine_h2b6s_6dpf/fish1/TSeries-titration-192trial-050"
        # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
        "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
    ](;resources...,
        registration_type = :dont_run
    );
##
# @pun (fish_dir, zseries) = recording1;
@pun (fish_dir, multimap_920) = recording1;
# fixedname = joinpath(fish_dir, "zseries.nrrd")
# save(fixedname, zseries)
movingname = joinpath(fish_dir, "multimap_920.nrrd")
save(movingname, multimap_920)
movingname