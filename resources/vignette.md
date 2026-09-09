# Workflow for: Targeted Compass Anaylsis of Pseudobulked Single-cell Datasets

This workflow enables Metabolic analysis of single-cell datasets.

Specifically, it uses the [Compass](https://github.com/wagnerlab-berkeley/Compass) algorithm to infer reaction activities of reactions related to Lipid, Amino Acid, and Central Carbon metabolism.
For computational efficiency, data are first psuedobulked in order to reduce complextiy from tens or hundreds of thousands of cells down to hundreds or thousands of pseudobulks.

## Before running this workflow

Before running this workflow, it is recommended to gain some familiarity with the dataset that you intend to analyze, and it is required that you obtain a license for a tool used inside the algorithm.

### Gaining familiarity with the single-cell data

One recommended way to become more familiar with the data you want to analyze, is to explore the data by analyzing it with the `scViz` workflow.
This workflow allows plotting of umaps, violins, dotplots, and compositional bar plots that can help build an understanding of the structures and makeup of cells' metadata, e.g. their annotations, identities, and other meaningful features.

To run that workflow, return to your Vulcan Project dashboard and create a workspace for that workflow.

### Obtaining your License to run this workflow

The proprietary Gurobi Optimizer is an inherent component that speeds up the Compass algorithm.  Luckily, Academic licenses are free.

Obtain your 'gurobi.lic' by following the instructions at [this page](https://support.gurobi.com/hc/en-us/articles/13232844297489-How-do-I-set-up-a-Web-License-Service-WLS-license).

You will provide the contents of this file to one of the early workflow inputs.

Your license should look something like this. (Note that the lines starting with "# " can be omitted.)
```
# Gurobi WLS license file
# Your credentials are private and should not be shared or copied to public repositories.
# Visit https://license.gurobi.com/manager/doc/overview for more information.
WLSACCESSID=########-####-####-####-############
WLSSECRET=########-####-####-####-############
LICENSEID=#######
```

## The Workflow in a bit more detail

The workflow starts with asking a few key parameterizations:
- dataset selection
- pseudobulking methodology
- some compass & post-processing parameters

After downloading the dataset into your workspace, you will then be able to select what metadata to use for pseudobulking.

### Module-Compass

The [Module-Compass](https://compass-wagnerlab.readthedocs.io/en/latest/module_compass.html) algorithm is integrated into the Vulcan workflow. Broadly speaking, Module-Compass is an algorithm that partitions the metabolic network into individual subsystems, then runs the Compass flux-balance analysis calculation on reactions within these subsystems. This provides orders-of-magnitude speedup compared to the original Compass algorithm by running calculations on smaller, discrete subsystems that preserve the topology of the original metabolic network. By default, Module-Compass supports the following subsystems:

- Central Carbon Metabolism
- Amino Acid Metabolism
- Lipid Metabolism

**Note that running the compass calculation takes a long time, often over a day!**


Lastly, you will be able to establish definitions for two groupings of pseudobulks to compare against each other.

An interprettable plot showing differential reaction scores between these groups will then be generated.

Raw compass outputs, as well as statistical outputs, can then be downloaded for local followup if desired.

### Normalization Choice

By default, reactions are normalized against counts of all genes. This preserves each cell's overall metabolic activity level, so that cells with globally higher biosynthetic activity (most notably proliferating cells) will carry this difference into the Compass analysis.


Normalizing against metabolic genes only removes this overall-activity signal, preserving relative usage across metabolic pathways. This is useful when the goal is to identify metabolic rewiring between groups independent of differences in overall metabolic activity.

### Reaction Metadata

To inspect the reactions in more detail, you can access the [Human1](https://github.com/wagnerlab-berkeley/Compass/tree/master/compass/Resources/Metabolic%20Models/Human1) or [Mouse1](https://github.com/wagnerlab-berkeley/Compass/tree/master/compass/Resources/Metabolic%20Models/Mouse1) resources directory of Compass where there are several .csv files that include metadata for the reactions. Alternatively, you can also visit [Metabolic Atlas](https://metabolicatlas.org/) to visualize the metabolic network.

### Deeper plot explanations

Each dot represents a reaction, colored red if higher in Group 1 (Cohen's d of 0 or above) or blue if higher in Group 2. Significant reactions (adjusted p below 0.05) are shown in solid color, while non-significant reactions are shown in lighter color. Reaction-level p-values are computed using the unpaired Wilcoxon rank-sum test (equivalent to the Mann-Whitney U test), comparing reaction consistency scores between the two groups. Effect sizes are calculated by Cohen's d, and p-values are FDR-corrected across all tested reactions.


Each triangle shows each subsystem's mean Cohen's d across its reactions. The triangle is colored red if its subsystem is significantly enriched for Group 1-higher reactions, blue if significantly enriched for Group 2-higher, and black if neither side is significant. Enrichment significance is computed with a hypergeometric test for enrichment of significant reactions within each subsystem, followed by FDR correction across subsystems.


To explore individual reactions further, please refer to the reaction-level statistics table within the downloadable reaction stats output. Please note that depending on the metabolic model, a reaction may appear in more than one pathway.