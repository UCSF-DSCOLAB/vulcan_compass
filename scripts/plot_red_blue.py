import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, hypergeom
from statsmodels.stats.multitest import multipletests
from matplotlib import __version__ as matplotlibversion
import json
if matplotlibversion < "3.4":
    print("Matplotlib versions older than 3.4 may not be able to generate figure 2E, as they do not support alpha arrays")

### Use Inputs
group_A_cells = list(pd.read_csv(snakemake.input['group_1_inds'], header=None)[0])[1:]
group_B_cells = list(pd.read_csv(snakemake.input['group_2_inds'], header=None)[0])[1:]
# Parse reaction file target
subsystem_full = snakemake.config['post_process_meta_subsystem']
subsystem = "carbon" if subsystem_full=="CENTRAL_CARBON_META_SUBSYSTEM" \
    else "lipid" if subsystem_full=="LIPID_META_SUBSYSTEM" \
    else "AA" if subsystem_full=="AA_META_SUBSYSTEM" \
    else "NOT FOUND ERROR"
reaction_suffix = "norm_sum" if snakemake.config['post_process_norm_method'] == "Sum__divide_by_sum_per_susbsystem" \
    else "norm_rank" if snakemake.config['post_process_norm_method'] == "Rank__rank_reactions_per_pseudobulk" \
    else "scores"
reactions_use = pd.read_csv(snakemake.input[f'{subsystem}_reaction_{reaction_suffix}'], sep="\t", index_col = 0)
reaction_metadata = pd.read_csv(f'output/compass_output/meta_subsystem_models/{subsystem_full}/{subsystem_full}_rxn_meta.csv', index_col = 0)

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

# Statistical test
# Unpaired Wilcoxon rank-sum test (equivalent ro Mann-Whitney U test)
# Positive values indicate higher potential activity in group_A cells than group_B cells.
wilcox_results = wilcoxon_test(reactions_use, group_A_cells, group_B_cells)
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

data.to_csv(snakemake.output['reaction_stats_csv'], sep='\t')

data_sig = data[data['adjusted_pval'] < 0.05].copy()
data_sig_pos = data[(data['adjusted_pval'] < 0.05) & (data['cohens_d'] > 0)].copy()
data_sig_neg = data[(data['adjusted_pval'] < 0.05) & (data['cohens_d'] < 0)].copy()

subsystems = data['SUBSYSTEM'].value_counts().index
hyper_pvals_pos = []
hyper_pvals_neg = []

for subsys in subsystems:
    #Positive
    M = data.shape[0]
    n = data_sig_pos.shape[0]
    N = data[data['SUBSYSTEM'] == subsys].shape[0]
    x = data_sig_pos[data_sig_pos['SUBSYSTEM'] == subsys].shape[0]
    pval_pos = hypergeom.sf(x-1, M, n, N)
    hyper_pvals_pos.append(pval_pos)
    #Negative
    M = data.shape[0]
    n = data_sig_neg.shape[0]
    N = data[data['SUBSYSTEM'] == subsys].shape[0]
    x = data_sig_neg[data_sig_neg['SUBSYSTEM'] == subsys].shape[0]
    pval_neg = hypergeom.sf(x-1, M, n, N)
    hyper_pvals_neg.append(pval_neg)

hyper_df_pos = pd.DataFrame({'pval': hyper_pvals_pos}, index=subsystems)
hyper_df_pos['adjusted_pval'] = multipletests(hyper_df_pos['pval'], method='fdr_bh')[1]

hyper_df_neg = pd.DataFrame({'pval': hyper_pvals_neg}, index=subsystems)
hyper_df_neg['adjusted_pval'] = multipletests(hyper_df_neg['pval'], method='fdr_bh')[1]

pd.DataFrame(
    {
        'Group1_higher_hypergeom_pval': hyper_pvals_pos,
        'Group1_higher_adjusted_pval': hyper_df_pos['adjusted_pval'],
        'Group2_higher_hypergeom_pval': hyper_pvals_neg,
        'Group2_higher_adjusted_pval': hyper_df_neg['adjusted_pval']
    },
    index=subsystems
).to_csv(snakemake.output['subsystem_stats_csv'], sep='\t')

### Plot
plt.figure(figsize=(12,12))
axs = plt.gca()
#Sorts the reactions for plotting
d = data.groupby('SUBSYSTEM')['cohens_d'].mean()
color = data['cohens_d'].map(lambda x: '#E31A1C' if x >= 0 else '#1F78B4')
size = np.array([100 if color[i] == 'orange' else 30 for i in range(len(color))])
alpha = data['adjusted_pval'].map(lambda x: 1.0 if x < 0.05 else 0.25)
d_color_pos = hyper_df_pos['adjusted_pval'].map(lambda x: '#FC9272' if x < 0.05 else 'black')
d_color_neg = hyper_df_neg['adjusted_pval'].map(lambda x: '#9ECAE1' if x < 0.05 else 'black')
d_color = pd.Series([
    'black' if s1 == 'black' and s2 == 'black' else (s1 if s2 == 'black' else s2)
    for s1, s2 in zip(d_color_pos, d_color_neg)
], index=d_color_pos.index)

sorted_subsystems = d[np.argsort(d)].index
axs.scatter(d[sorted_subsystems], d[sorted_subsystems].index, marker='^', c=d_color[sorted_subsystems], s=50)
axs.scatter(data[data['SUBSYSTEM'].isin(sorted_subsystems)]['cohens_d'].values,
            data[data['SUBSYSTEM'].isin(sorted_subsystems)]['SUBSYSTEM'].values,
            c=color, alpha=alpha, s=size)
axs.scatter(d[sorted_subsystems], d[sorted_subsystems].index, marker='^', c=d_color[sorted_subsystems], s=90, edgecolors='black', label='mean of all reactions')

axs.set_xlabel(f"(higher in Group 2)                          Cohen's d                          (higher in Group 1)")

plt.legend()

plt.savefig(snakemake.output['plot'], dpi=300, bbox_inches="tight")
