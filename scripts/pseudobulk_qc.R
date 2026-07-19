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

df <- df[,c(
    "metabolism_counts_fraction",
    "avg_nCounts_per_cell__all_genes",
    "avg_nCounts_per_cell__metab_targets",
    ct_col
)]
names_use <- c(
    "Metabolism Targets Fraction",
    "Total Counts - All Genes",
    "Total Counts - Metabolic Target Genes"
)
colnames(df)[1:3] <- names_use 

# Median ranges
df_lines <- data.frame(
    var.which = names_use
)
df_lines$median <- sapply(
    names_use, function(col) {
        median(df[,col])
    }
)
df_lines$IQR <- sapply(
    names_use, function(col) {
        IQR(df[,col])
    }
)


cts <- unique(df[,ct_col])
n_cts <- length(cts)
n_char_max_ct <- max(nchar(cts))

### Plot!
ts_log('Plotting')
width <- 0.7 + 0.3*n_cts
height <- 2.5*3 + 0.075*n_char_max_ct
png(output_path('all_plots'), w = width*75, h = height*75, res=75)
yPlot(
    df,
    names_use,
    ct_col,
    main = "Metabolism Count Metrics for Normalization Choice Assessment",
    sub = "Data points represent each reatined pseudobulk\nValues represent the mean across pseudobulks' contituent cells",
    split.ncol = 1,
    split.adjust = list(scale = 'free_y'),
    plots = c("vlnplot", "jitter"),
    vlnplot.quantiles = c(0.25, 0.5, 0.75),
    vlnplot.lineweight = 0.5,
    vlnplot.scaling = "width",
    jitter.width = 0.8,
    legend.show = FALSE) +
    geom_hline(
        data = df_lines,
        mapping = aes(y=median),
        linetype = 'dashed'
    ) +
    geom_hline(
        data = df_lines,
        mapping = aes(y=median - 0.75 * IQR),
        lineweight = 0.5,
        linetype = 'dotted'
    ) +
    geom_hline(
        data = df_lines,
        mapping = aes(y=median + 0.75 * IQR),
        lineweight = 0.5,
        linetype = 'dotted'
    )
dev.off()

ts_log('DONE')
