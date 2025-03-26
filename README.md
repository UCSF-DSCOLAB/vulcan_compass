# vulcan_compass
This repo encompasses a Snakemake workflow and additional vulcan configuration files that allow running the [Compass pipeline](https://github.com/wagnerlab-berkeley/Compass) on single-cell datasets of the UCSF Data Library.

The workflow is currently in development, but the plan is to enable:
1. Pseudobulking data by sample and cell type
2. Running Compass targeting core central carbon metabolism, amino accid metabolism, and lipid metabolism
   1. Choosing between multiple options for data normalization
3. Running statistical analysis on Compass results between groups, and visualizing this output
4. Running DGE, with DESeq2, then plotting volcanos with metabolism reaction-related genes labeled for validation 
