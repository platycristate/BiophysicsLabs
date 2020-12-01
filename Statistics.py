# threshold p-value for t_test
thresh = 0.05


# Calculates significance of difference between groups
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import shapiro, ttest_rel
plt.style.use('dark_background')

experiment_name = 'Task2'

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
plt.savefig(experiment_name + '/Stats/histograms_amps.pdf')

hist = areas.hist(bins=15, figsize=(12,10));
plt.savefig(experiment_name + '/Stats/histograms_areas.pdf')

hist = hws.hist(bins=15, figsize=(12,10));
plt.savefig(experiment_name + '/Stats/histograms_hws.pdf')


# t-test 
keys = list(amplitudes_sep.keys())
vals = amplitudes_sep.values()
p_values = []
mv = min([len(v) for v in vals])
for i in range(len(keys)):
    for j in range(i+1, len(keys)):
        p_values.append( (str(keys[i]), str(keys[j]), ttest_rel(amplitudes_sep[keys[i]][:mv], amplitudes_sep[keys[j]][:mv])[1]))

pval_df = pd.DataFrame(columns=keys, index=keys)
for i in p_values:
	if i[2] < thresh:
		pval_df.loc[i[0], i[1]] = '<' + str(thresh)
pval_df.to_csv(experiment_name + '/Stats/t_tests_pvals.csv')
