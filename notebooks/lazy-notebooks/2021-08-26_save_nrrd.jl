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
recording = Recordings[
        "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    ](;resources...,
        registration_type = :dont_run
    );
##
@pun zbrain_registered = recording