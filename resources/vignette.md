# Workflow for: Targeted Compass Anaylsis of Pseudobulked Single-cell Datasets

This workflow enables differential analysis of Lipid, Amino Acid, and Central Carbon Metabolic Reaction Activities across samples after psuedobulking a single-cell RNAseq dataset.

## Obtaining your License to run this workflow

The proprietary Gurobi Optimizer is an inherent component that speeds up the Compass algorithm.  Luckily, Academic licenses are free.

Obtain your 'gurobi.lic' by following the instructions at [this page](https://support.gurobi.com/hc/en-us/articles/13232844297489-How-do-I-set-up-a-Web-License-Service-WLS-license). Then provide the contents of this file to a workflow input.

## The Workflow in a bit more detail

The workflow starts with asking a few key parameterizations:
- dataset selection
- pseudobulking methodology
- some compass & post-processing parameters

After downloading the dataset into your workspace, you will then be able to select what metadata to use for pseudobulking.

It then performs pseudobulking and runs Compass towards processes in the categories:
- Central Carbon Metabolism
- Amino Acid Metabolism
- Lipid Metabolism

**Note that running the compass calculation takes a long time, often over a day!**

Lastly, you will be able to establish definitions for two groupings of pseudobulks to compare against eachother,
and an interprettable plot will be generated.

Raw compass outputs as well as statistical outputs can then be downloaded for local followup if desired.

