suppressPackageStartupMessages({
    library(dittoViz)
    library(dataflow)
})

### Function imports from essential scripts
ts_log <- function(..., sep = "") {
    cat("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), "] ", paste(..., sep = sep), "\n", sep = "")
}

### Parameters
ts_log('Reading in data and inputs')
norm_method <- input_str('pre_process_norm_method')
samp_col <- input_str('sample_id_column')
ct_col <- input_str('cell_type_column')
df <- read.table(input_path('pseudo_metadata'), sep = "\t", header = TRUE, row.names = 1)

ts_log('Reading in data and inputs')
cts <- unique(df[,ct_col])
n_cts <- length(cts)
n_char_max_ct <- max(nchar(cts))
med_n_samp <- median(table(df[,samp_col], df[,ct_col]))

### Establish some plot features
ts_log('Prepping to plot')
reps <- ifelse(
    med_n_samp <= 20,
    c("jitter", "boxplot"),
    c("jitter", "vlnplot")
)
width <- 0.7 + 0.3*n_cts
height <- 2.5*3 + 0.075*n_char_max_ct

### Plot!
ts_log('Plotting')
png(output_path('all_plots'), w = width*75, h = height*75, res=75)
yPlot(
    df,
    c('metabolism_counts_fraction', 'avg_nCounts_per_cell__all_genes', 'avg_nCounts_per_cell__metab_targets'),
    ct_col,
    main = "Metabolism Ammount Comparison Metrics",
    sub = "per sample/pseudobulk, grouped by cell type",
    split.ncol = 1,
    split.adjust = list(scale = 'free_y'),
    plots = reps,
    vlnplot.quantiles = c(0.25, 0.5, 0.75),
    vlnplot.lineweight = 0.5,
    vlnplot.scaling = "width",
    boxplot.width = 0.8,
    boxplot.fill = FALSE,
    boxplot.lineweight = 0.5,
    jitter.width = 0.8,
    legend.show = FALSE)
dev.off()

ts_log('DONE')
