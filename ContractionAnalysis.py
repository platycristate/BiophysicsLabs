#--------------------------------------INPUT INFORMATION-------------------------
###################################################################################


experiment_name = 'Task2'
path_to_file = 'task2.csv'


substances = {
            'Cheler': (43, 53),
            'Oxyt': (84, 92),
            'Oxyt+Cheler': (124, 134)
            }

Δt=0.1/60

peaks_distance = 200
width_height = 0.93

#################################################################################


#------------------------------------- GENERAL SETTINGS--------------------------
vals = list(substances.values())
vals.sort(key=lambda x: x[0])
substances['control'] = (0, vals[0][0])
print(substances)
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
from BaselineRemoval import BaselineRemoval
plt.style.use('dark_background')

color=iter(cycle(plt.cm.rainbow(np.linspace(0,1,len(substances)))))

# loading the data
contrac = np.loadtxt(path_to_file, delimiter=',',
                    usecols=1, skiprows=12000)

# Determine duration
time = np.arange(0, Δt*len(contrac), Δt)

# Find baseline and remove it from contraction trace
base = peakutils.baseline(contrac, 5)
contrac_aligned = contrac - base

# General view after alignment
fig, axs = plt.subplots(2, 1, figsize=(20, 8), dpi=130)

axs[0].plot(time, contrac)
axs[0].plot(time, base, color='red')
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

# Finding variables required for calculating amplitude, areas under peaks, half-widths
aligned = BaselineRemoval(contrac_aligned)
zhang = aligned.ZhangFit()

fig, axs =  plt.subplots(figsize=(20, 4), dpi=130)
axs.plot(time, zhang)
fig.savefig(experiment_name + '/ZhangFit.pdf')

# Determine threshold height of the peak
max_val = np.where(zhang == zhang.max())[0][0]
peaks_height = (zhang.max() - zhang[max_val-peaks_distance: max_val+peaks_distance].min()) * 0.2

peaks, peak_values = find_peaks(zhang, height=peaks_height,
        distance=peaks_distance)
widths, widths_heights, left, right = peak_widths(zhang, peaks,
        rel_height=width_height)
left = left.astype(int)
right = right.astype(int)


# to remove outliers; so called "temporary solution to the problem"
# Calculaitng half-widths 
half_widths, widths_heights_half, left_half, right_half = peak_widths(zhang, peaks, rel_height=0.5)
half_widths *= Δt*60000

fig, axs = plt.subplots(len(substances)+1, 1, figsize=(20, 4*(len(substances)+1)), dpi=130)
axs[0].plot(time, contrac_aligned)
axs[0].plot(time[peaks], contrac_aligned[peaks], 'x')

for ax, substance, times in zip(axs[1:], list(substances.keys()), list(substances.values())):
	ax.plot(time, contrac_aligned, label=substance)
	ax.legend()
	ax.hlines(contrac_aligned[left_half.astype(int)], left_half*Δt, right_half*Δt, color='red')
	ax.set_xlim(list(times))

axs[-1].set(xlabel='time [min]')
fig.savefig(experiment_name + '/Peaks_half_widths_amplitude.pdf')

# plotting limits of integratiobn
fig, axs = plt.subplots(len(substances), 1,
						figsize=(20, 4*len(substances)), dpi=130);

for ax, substance, times in zip(axs, list(substances.keys()), list(substances.values())):
	ax.plot(time, contrac_aligned, label=substance)
	ax.legend()
	ax.hlines(widths_heights, left*Δt, right*Δt, color='red')
	ax.set_xlim(list(times))

axs[-1].set(xlabel='time [min]')

fig.savefig(experiment_name + '/Limits_of_integration.pdf')

# find areas and amplitudes
amps = []
areas = []
for idx, peak in enumerate(peaks):
    amp = contrac_aligned[peak]
    amps.append(amp)
    area = simps(
            contrac_aligned[left[idx]:right[idx]],
            60000*time[left[idx]:right[idx]])
    areas.append(area)

# Separating amplitudes, areas, half_widths
amplitudes_sep = {key: [] for key in substances.keys()}
areas_sep = {key: [] for key in substances.keys()}
half_widths_sep = {key: [] for key in substances.keys()}

for key, times in substances.items():
	amplitudes_sep[key] = amps[np.argmax(time[peaks] > times[0]) : np.argmax(time[peaks] > times[1])]
	areas_sep[key] = areas[np.argmax(time[peaks] > times[0]) : np.argmax(time[peaks] > times[1])]
	half_widths_sep[key] = half_widths[np.argmax(time[peaks] > times[0]) : np.argmax(time[peaks] > times[1])]

#------------------------------ PLOTTING ERRORBARS----------------------------------
# calculating means and standart errors
labels = [key for key in substances.keys()]
x = np.arange(len(labels))

for d in [amplitudes_sep, areas_sep, half_widths_sep]:
	for k, v in d.items():
		d[k] = (np.mean(v), np.std(v)/np.sqrt(len(v)))

fig, axs = plt.subplots(1, 3, dpi=150, figsize=(16, 4))
for idx, d in enumerate([amplitudes_sep, areas_sep, half_widths_sep]):
	y = [ d[ labels[i] ][0] for i in range(len(labels))]
	error = [ d[ labels[i] ][1] for i in range(len(labels))]
	axs[idx].bar(x, y, yerr=error, align='center', alpha=0.8, ecolor='red', capsize=10, facecolor='silver')
	axs[idx].set_xticks(x)
	axs[idx].set_xticklabels(labels)

axs[0].set_title('Contraction amplitude, g')
axs[2].set_title('Contraction half_widths, ms')
axs[1].set_title('Contraction areas, $g \cdot ms$')

fig.savefig(experiment_name + '/Errorbars.pdf')





