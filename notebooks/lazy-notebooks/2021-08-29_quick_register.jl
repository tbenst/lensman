using Lensman
import Lensman: @pun
resource = Resources()
r = Recordings[
    # "2021-07-14_rsChRmine_h2b6s_5dpf/fish2/TSeries-lrhab-118trial-069"
    "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"
    ](
    ;resource...);
##
@pun (multimap_ants_cmd, zseries, multimap_820, oir_dir) = r

ants_register(zseries, multimap_820[:,:,3,1,:]; interpolation = "WelchWindowedSinc",
            histmatch = 0, sampling_frac = 0.25, maxiter = 200, threshold=1e-8,
            use_syn = true, synThreshold = 1e-7, synMaxIter = 200,
            save_dir=oir_dir, run = false)
