using Lensman
resources = Resources()
recording1 = Recordings[
        # "2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-titration-192trial-062"
        "2021-06-08_rsChRmine_h2b6s/fish2/TSeries-lrhab-titration-123"
    ](;resources...,
    );
##
trial_order = recording1[:trial_order];
stim_start_idx = recording1[:stim_start_idx];
stim_end_idx = recording1[:stim_end_idx];
@show (trial_order, stim_start_idx, stim_end_idx)