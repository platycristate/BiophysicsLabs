#----------------------------------------INPUT INFORMATION------------------------------

# threshold p-value for t_test or U-test
thresh = 0.05
experiment_name = 'Task2'


#---------------------------------------------------------------------------------------

import pickle
import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
from scipy.stats import shapiro, ttest_rel, mannwhitneyu
plt.style.use('dark_background')


# loading the data
with open(experiment_name + '/amplitudes_sep.p', 'rb') as f:
    amplitudes_sep = pickle.load(f)

with open(experiment_name + '/areas_sep.p', 'rb') as f:
    areas_sep = pickle.load(f)

with open(experiment_name + '/half_widths_sep.p', 'rb') as f:
    half_widths_sep = pickle.load(f)

# histograms
amps = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in amplitudes_sep.items() ]))
areas = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in areas_sep.items() ]))
hws = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in half_widths_sep.items() ]))


hist = amps.hist(bins=15, figsize=(12,10));
plt.suptitle('Histograms for amplitudes', fontsize=20)
plt.savefig(experiment_name + '/Stats/histograms_amps.pdf')

hist = areas.hist(bins=15, figsize=(12,10));
plt.suptitle('Histograms for areas', fontsize=20)
plt.savefig(experiment_name + '/Stats/histograms_areas.pdf')

hist = hws.hist(bins=15, figsize=(12,10));
plt.suptitle('Histograms for half-widths', fontsize=20)
plt.savefig(experiment_name + '/Stats/histograms_hws.pdf')


# t-test 
keys = list(amplitudes_sep.keys())
vals = amplitudes_sep.values()
p_values = []

# t-test requires the same amount of points
# for each sample =(
mv = min([len(v) for v in vals])
# tests for normal distribution
# make random samples of length of the smallest group
amplitudes_sep = {key: random.sample(val, mv) for key, val in amplitudes_sep.items()}
shapiro_res = {key: shapiro(val)[1] for key, val in amplitudes_sep.items()}

for i in range(len(keys)):
    for j in range(i+1, len(keys)):
        if shapiro_res[keys[i]] > thresh and shapiro_res[keys[j]] > thresh:
            p_values.append( (str(keys[i]), str(keys[j]), ttest_rel(amplitudes_sep[keys[i]], amplitudes_sep[keys[j]])[1]) )
        else:
            # Mann-Whitney rank test
            p_values.append( (str(keys[i]), str(keys[j]), mannwhitneyu(amplitudes_sep[keys[i]], amplitudes_sep[keys[j]])[1]) )

# Bonferoni correction	
thresh /= len(p_values) # divide by the number of comparisons

# Find pairs with significant difference in their means
pval_df = pd.DataFrame(columns=keys, index=keys)
for i in p_values:
	if i[2] < thresh:
		pval_df.loc[i[0], i[1]] = '<' + str(round(thresh, 4))
pval_df.to_csv(experiment_name + '/Stats/t_tests_pvals.csv')

# saving results of shapiro test
print('Results of Shapiro-Wilco tests: \n', shapiro_res)

with open(experiment_name + '/Stats/shapiro_res.p', 'wb') as f:
	pickle.dump(shapiro_res, f)
