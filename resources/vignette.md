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

It then performs pseudobulking and runs Compass towards processes in the categories:
- Central Carbon Metabolism
- Amino Acid Metabolism
- Lipid Metabolism

**Note that running the compass calculation takes a long time, often over a day!**

Lastly, you will be able to establish definitions for two groupings of pseudobulks to compare against eachother.
An interprettable plot showing differential reaction scores between these groups will then be generated.

Raw compass outputs, as well as statistical outputs, can then be downloaded for local followup if desired.

