experiment_name = '24_November_2020'
path_to_file = 'contraction.csv'
duration = 110 # mins
substances = {
			'KCL': (30, 37), 
			'Carbachol': (56, 62)
			 }

# Create directory for storing the results
import os
if not os.path.exists(experiment_name):
	os.mkdir(experiment_name)

	
#---------------------------- ANALYSIS PART-----------------------------------
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import peakutils
from scipy.signal import find_peaks, peak_widths
from scipy.integrate import simps
plt.style.use('dark_background')

color=iter(cycle(plt.cm.rainbow(np.linspace(0,1,len(substances)))))

# loading the data
contrac = np.loadtxt(path_to_file, delimiter=',')

# Determine time vector
Δt = duration / len(contrac)
time = np.arange(0, duration, Δt)

# Find baseline and remove it from contraction trace
base = peakutils.baseline(contrac, 5)
contrac_aligned = contrac - base

# General view after alignment
fig, axs = plt.subplots(2, 1, figsize=(20, 8), dpi=130)

axs[0].plot(time, contrac)
for substance, times in substances.items():
	axs[0].axvspan(*times, facecolor=next(color), edgecolor='None', alpha=0.5, label=substance)
	
axs[0].legend(title='Application')
axs[0].set_title('Before the alignment procedure')

axs[1].plot(time, contrac_aligned)
for substance, times in substances.items():
	axs[1].axvspan(*times, facecolor=next(color), edgecolor='None', alpha=0.5, label=substance)

axs[1].legend(title='Application')
axs[1].set_title('After the alignment procedure')

axs[1].set(xlabel = 'time [min]')

fig.savefig(experiment_name + '/General_traces.pdf')
			
			
			
			