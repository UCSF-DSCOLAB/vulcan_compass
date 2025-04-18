import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from matplotlib import __version__ as matplotlibversion
if matplotlibversion < "3.4":
    print("Matplotlib versions older than 3.4 may not be able to generate figure 2E, as they do not support alpha arrays")

import os
os.chdir('/krummellab/data1/danb/compass')

### Tweakable Parameters
subsystem_target = input['subsystem']
group_defs = input['group_defs']

### Standard parameters
reactions_file = f'output/compass_outputs/{subsystem_target}/reactions.tsv'
reaction_metadata_file = f'output/compass_outputs/{subsystem_target}/{subsystem_target}_rxn_meta.csv'
sample_metadata_file = input['pseudo_metadata']

### Functions
def cohens_d(x, y):
    pooled_std = np.sqrt(((len(x)-1) * np.var(x, ddof=1) 
                          + (len(y)-1) * np.var(y, ddof=1)) / 
                             (len(x) + len(y) - 2))
    return (np.mean(x) - np.mean(y)) / pooled_std
def wilcoxon_test(consistencies_matrix, group_A_cells, group_B_cells):
	"""
		Performs an unpaired wilcoxon test (or mann-whitney U test) for each reaction between group_A and group_B
	"""
	#per reaction/meta-reaction, conduct wilcoxon test between group_A and group_B
	group_A = consistencies_matrix.loc[:,group_A_cells]
	group_B = consistencies_matrix.loc[:,group_B_cells]
	results = pd.DataFrame(index = consistencies_matrix.index, columns = ['wilcox_stat', 'wilcox_pval', 'cohens_d'], dtype='float64')
	for rxn in consistencies_matrix.index:
		A, B = group_A.loc[rxn].to_numpy().ravel(), group_B.loc[rxn].to_numpy().ravel()
		#sometimes there's a solitary value, and we don't want to test then
		if len(np.unique(A)) == 1 and len(np.unique(B)) == 1:
			if np.unique(A) == np.unique(B):
				#we've got no data. set p-value to 1 and skip!
				#(p-value needs to be 1 so multipletests doesn't cry)
				results.loc[rxn, ['wilcox_pval']] = 1
				continue
		stat, pval = mannwhitneyu(A, B, alternative='two-sided')
		c_d = cohens_d(A, B)
		results.loc[rxn, ['wilcox_stat', 'wilcox_pval', 'cohens_d']] = stat, pval, c_d
	results['adjusted_pval'] = np.array(multipletests(results['wilcox_pval'], method='fdr_bh')[1], dtype='float64')
	return results
def get_reaction_consistencies(compass_reaction_penalties, min_range=1e-3):
    """
        Converts the raw penalties outputs of compass into scores per reactions where higher numbers indicate more activity
    """
    df = -np.log(compass_reaction_penalties + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
    df = df - df.min().min()
    return df
def reaction_out_index_to_id(index, reaction_metadata):
    # Goal: Remove '_pos' or '_neg' and ensure matches with known reactions
    if index in reaction_metadata.index:
        return index
    # without last 4 chars, a.k.a. the '_pos' or '_neg' (all)
    elif r[:-4] in reaction_metadata.index:
        return index[:-4]
    else:
        raise Exception("reaction unknown")

##### Primary Processing

### Loading needed data files
reaction_penalties = pd.read_csv(reactions_file, sep="\t", index_col = 0)
cell_metadata = pd.read_csv(sample_metadata_file, index_col = 0)
reaction_metadata = pd.read_csv(reaction_metadata_file, index_col = 0)

### --- Quirk for this data --- ###
# Oops... re-making from the pseudobulk's names
import re
def empty_or_match(search):
    if search is None:
        return ''
    else:
        return search.group(1)
cell_metadata['manual.celltype'] = [empty_or_match(re.search(r"^.*__(.+)$", i)) for i in cell_metadata.index]
### --- End Quirk --- ###

# Gather samples in comparator groups
################################################
# TBD From here on... need to implement the reciprocal of js group constraint selection UI to parse groupings
################################################
group_A_cells = cell_metadata.index[cell_metadata[diff_btwn_column] == diff_group_A]
group_B_cells = cell_metadata.index[cell_metadata[diff_btwn_column] == diff_group_B]

# Transform reaction penalties per cell/sample into scores that are higher the more active the reaction is predicted to be
reaction_consistencies = get_reaction_consistencies(reaction_penalties)

# Statistical test
# Unpaired Wilcoxon rank-sum test (equivalent ro MAnn-Whitney U test)
# Positive values indicate higher potential activity in group_A cells than group_B cells.
wilcox_results = wilcoxon_test(reaction_consistencies, group_A_cells, group_B_cells)
wilcox_results['metadata_r_id'] = ""
for r in wilcox_results.index:
    wilcox_results.loc[r, 'metadata_r_id'] = reaction_out_index_to_id(r, reaction_metadata)

# Merge wilcox_results with reaction info
W = wilcox_results.merge(reaction_metadata, how='left', 
                         left_on='metadata_r_id', right_index=True, validate='m:1')
# Apply some filtering
W = W[~W['EC-NUMBER'].isna()]
W.loc[(W['EQUATION'].map(lambda x: '[m]' not in x)) & (W['SUBSYSTEM'] == "Citric acid cycle"), 'SUBSYSTEM'] = 'Other'

## Final data prep
# Trim reactions from certain subsystems
data = W[~W['SUBSYSTEM'].isin(["Miscellaneous", "Unassigned", "Other"])]
data = data[~data['SUBSYSTEM'].map(lambda x: "Transport" in x or "Exchange" in x)]
# Trim subsystems with fewer than 5 measurements
items, counts = np.unique(data['SUBSYSTEM'], return_counts=True)
items = [items[i] for i in range(len(items)) if counts[i] > 5] #filter(n() > 5) %>%
data = data[data['SUBSYSTEM'].isin(items)]

### Plot
plt.figure(figsize=(12,12))
axs = plt.gca()
#Sorts the reactions for plotting
d = data[data['adjusted_pval'] < 0.1].groupby('SUBSYSTEM')['cohens_d'].median().abs()
axs.scatter(d[d.argsort], d[d.argsort].index, alpha=0)
color = data['cohens_d'].map(lambda x: 'r' if x >= 0 else 'b')
alpha = data['adjusted_pval'].map(lambda x: 1.0 if x < 0.1 else 0.25)
axs.scatter(data['cohens_d'], data['SUBSYSTEM'], c=color, alpha=alpha)
axs.set_xlabel("Cohen's d")
plt.savefig(output['plot'], dpi=300, bbox_inches="tight")
