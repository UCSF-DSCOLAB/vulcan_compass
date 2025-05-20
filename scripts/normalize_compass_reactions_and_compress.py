from archimedes.functions.dataflow import input_path, output_var, output_tsv, output_tgz
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from snakemake.script import snakemake


### First: Normalizations
def get_reaction_consistencies(compass_reaction_penalties: pd.DataFrame, min_range=1e-3):
    """
        Converts the raw penalties outputs of compass into scores per reactions where higher numbers indicate more activity
    """
    df = -np.log(compass_reaction_penalties + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
    df = df - df.min().min()
    return df

readme = "Files in this directory:\n\n- reactions.tsv: Raw Penalties output from compass where lower values equate to higher activity\n- reaction_scores.tsv: Transformation and minor trimming of raw penalties into \"Scores\" where higher numbers indicate higher activity. Steps: 1. Transformation: -log(x+1). 2. Trimming: Reactions removed if trasformed range, across all samples, is less than or equal to 0.001. 3. Final Transformation: x - min(all_x) to have at least one 0 value in the tsv.\n- reactions_norm_sum.tsv: Normalization by dividing each reaction score by the per-column (per-sample) sum of reaction scores.\n- reactions_norm_rank.tsv: Noramlization by ranking reaction scores per column, where the lowest value becomes 1, next lowest becomes 2, and so on.\n\nNote: Normalization of reaction scores is not absolutely required, especially when the expression data was trimmed to only metabolic genes and then renormalized before compass was run. That said, the rank based normalization seems to be the prefered normalization method due to magnitudinal differences in ranges of reactions' scores."

def norm_subsystem(in_: str, out_score: str, out_sum: str, out_rank: str, out_readme: str):
    df_in = pd.read_csv(input_path(in_, snakemake), sep="\t", index_col = 0)
    df_score = get_reaction_consistencies(df_in)
    output_tsv(df_score, out_score, snakemake)
    df_sum = df_score / df_score.sum(axis=0)
    output_tsv(df_sum, out_sum, snakemake)
    df_rank = np.apply_along_axis(rankdata, axis=0, arr=df_score)
    output_tsv(df_rank, out_rank, snakemake)
    output_var(readme, out_readme, snakemake)

norm_subsystem("carbon_reactions", "carbon_reaction_scores", "carbon_reaction_norm_sum", "carbon_reaction_norm_rank", "carbon_readme")
norm_subsystem("lipid_reactions", "lipid_reaction_scores", "lipid_reaction_norm_sum", "lipid_reaction_norm_rank", "lipid_readme")
norm_subsystem("AA_reactions", "AA_reaction_scores", "AA_reaction_norm_sum", "AA_reaction_norm_rank", "AA_readme")

### Last: compress the folder
output_tgz("output/compass_output", "compass_tgz", snakemake)



