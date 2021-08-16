using Lensman, Images, Glob, NPZ, PyCall, DataFrames, ImageSegmentation, 
      Random, Gadfly, Statistics, PyCall, ProgressMeter, HDF5, Distributed
import PyPlot
import Lensman: @pun, @assign
plt = PyPlot
matplotlib = plt.matplotlib
np = pyimport("numpy")
import Plots
import Plots: heatmap
using Thunks
import Base.Threads: @spawn, @sync, @threads
L = Lensman
init_workers()
##
r = Recordings["2021-07-14_rsChRmine_h2b6s_5dpf/fish1/TSeries-lrhab-118trial-061"]();
@pun (tseries, tseriesT) = r
##
@pun window_secs = r.settings
@pun (trial_average, vol_rate, window_len, recording_folder, fish_name,
    expName, tseriesZ, tseriesW, cells, exp_date, plot_dir, nStimuli,
    window_len, etl_vals, mm820_registered, mm920_registered
) = r
## idea: plot violin plot of df/f voxels in medulla vs optic tectum
