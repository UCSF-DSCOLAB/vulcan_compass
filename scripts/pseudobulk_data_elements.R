suppressPackageStartupMessages({
    library(Seurat)
    library(dataflow)
    library(edgeR)
})

### Parameters
scd_rds <- input_path('scd_rds')
min_cells <- input_num('min_cells')
samp_col <- input_str('sample_id_column')
ct_col <- input_str('cell_type_column')
norm_method <- input_str('pre_process_norm_method')
norm_genes <- input_str('target_genes')

### Function imports from essential scripts
logger_ts <- function(...) {
    cat("[ ", format(Sys.time()), " ] ", ..., "\n", sep = "")
}
logger <- function(...) {
    cat(..., "\n", sep = "")
}

### Data Input
logger_ts("Loading Data and Prepping")
sobj <- readRDS(scd_rds)
logger_ts('Done')

# Metadata, per-cell
meta <- sobj@meta.data

# Counts - delog
logger_ts("Delogging")
counts <- GetAssayData(sobj, assay = 'RNA', layer = 'counts')
rm(sobj)
gc(verbose=FALSE)

### Function for metadata collapse
meta_collapse <- function(df_col, name) {
    if (is.numeric(df_col)) {
        mean(df_col)
    } else {
        if (name %in% meta_ignored) {
            NA
        } else {
            if (length(unique(df_col))>1) {
                logger("\tadding ", name, " to ignored metadata columns")
                meta_ignored <<- c(meta_ignored, name)
                NA
            } else {
                df_col[1]
            }
        }
    }
}
### Function for summarizing metadata for later UIs
summarize_discrete_metadata <- function(metadata) {
    
    output <- list()
    for (i in colnames(metadata)) {
        this_data <- metadata[,i, drop = TRUE]
        if (!is.numeric(this_data)) {
            # Remove any empty levels, but keep order if already a factor
            used <- levels(as.factor(as.character(this_data)))
            levs <- levels(as.factor(this_data))
            this_out <- levs[levs %in% used]

            # 'key' as numbers to match with the pandas output format already being used
            names(this_out) <- seq_along(this_out)

            output[[i]] <- this_out
        }
    }

    output
}

### Pseudobulk Counts
# Cells
cell_annots <- meta[,ct_col]
celltypes <- unique(cell_annots)
celltypes <- celltypes[!is.na(celltypes)]
# Donors/Samples
cell_samples <- meta[,samp_col]
samples <- unique(cell_samples)
samples <- samples[!is.na(samples)]
# Tracking skippage due to min_cells
skipped_sample <- NULL
skipped_cell <- NULL
skipped_i <- NULL
i <- 0
# Tracking metadata skippage
meta_ignored <- NULL

logger_ts("Starting Pseudobulking")
pseudo_mat <- NULL
pseudo_meta <- NULL
for (ct in celltypes) {
    logger_ts("\tworking on ", ct)
    is_cell <- cell_annots==ct & !is.na(cell_annots)
    for (samp in samples) {
        # logger_ts("\t\tworking on ", samp)
        i <- i + 1
        is_set <- is_cell & cell_samples==samp & !is.na(cell_samples)
        if (sum(is_set) > min_cells) {
            # Trim matrix and make dense
            counts_set <- counts[,is_set]
            this_pseudo_mat <- rowSums(counts)
            meta_set <- meta[is_set,]
            this_pseudo_meta <- meta_set[1,,drop=FALSE]
            for (meta_i in 1:ncol(meta_set)) {
                this_pseudo_meta[1,meta_i] <- meta_collapse(meta_set[,meta_i], names(meta_set)[meta_i])
            }
            this_pseudo_meta$pseudobulk_cell_num <- sum(is_set)
            this_name <- paste0(samp, "__", ct)
            # Add column
            if (identical(pseudo_mat, NULL)) {
                pseudo_mat <- data.frame(this_pseudo_mat)
                colnames(pseudo_mat) <- this_name
                pseudo_meta <- meta[1,]
                pseudo_meta$pseudobulk_cell_num <- NA
                pseudo_meta[1,] <- this_pseudo_meta
                rownames(pseudo_meta) <- this_name
            } else {
                pseudo_mat <- cbind(pseudo_mat, this_pseudo_mat)
                colnames(pseudo_mat)[ncol(pseudo_mat)] <- this_name
                if (ncol(this_pseudo_meta) != ncol(pseudo_meta)) {
                    pseudo_meta <- pseudo_meta[,colnames(this_pseudo_meta), drop = FALSE]
                }
                pseudo_meta <- rbind(pseudo_meta, this_pseudo_meta)
                rownames(pseudo_meta)[nrow(pseudo_meta)] <- this_name
            }
        } else {
            skipped_sample <- c(skipped_sample, samp)
            skipped_cell <- c(skipped_cell, ct)
            skipped_i <- c(skipped_i, i)
            logger_ts("\t\tskipping ", samp, ", only ", sum(is_set), " cells")
        }
    }
}
rownames(pseudo_mat) <- rownames(counts)

logger_ts("Pseudobulking Complete")
if (length(skipped_i)>0) {
    logger("NOTE: ", length(skipped_i), " pseudobulks were skipped due to minimum cell threshold of ", min_cells, ":")
    logger("Samples with any skips, and number of celltypes of each that were skipped:")
    print(sort(table(skipped_sample)))
    logger("Celltypes with any skips, and number of celltypes of each that were skipped:")
    print(sort(table(skipped_cell)))
}
if (length(meta_ignored)>0) {
    logger("NOTE: ", length(meta_ignored), " discrete metadata columns were skipped in collapsing due to inconsistent values across pseudobulked cells:")
    logger("Metadata columns with inconsistencies:")
    logger("\t", paste0(meta_ignored, collapse = "\n\t"))
    pseudo_meta <- pseudo_meta[,!colnames(pseudo_meta)%in% meta_ignored]
}

logger_ts("Calculating scaling factors and metabolic gene percentages")
dgelist_all <- calcNormFactors(DGEList(counts = pseudo_mat), method = "TMM")

genes_in <- intersect(rownames(pseudo_mat), norm_genes)
metab_mat <- pseudo_mat[genes_in,]
dgelist_metab <- calcNormFactors(DGEList(counts = metab_mat), method = "TMM")

pseudo_meta$metabolism_counts_fraction <- colSums(metab_mat)/colSums(pseudo_mat)
pseudo_meta$scaling_factors__all_genes <- dgelist_all$samples$norm.factors
pseudo_meta$scaling_factors__metab_targets <- dgelist_metab$samples$norm.factors

if (norm_method=="Yes__Only_to_metabolic_genes") {
    logger_ts("Note: Focusing toward metabolic genes only")
    pseudo_norm <- cpm(dgelist_all, log = FALSE)
} else if (norm_method!="No__To_all_genes") {
    logger_ts("norm_method: ", norm_method)
    stop("WORKFLOW ERROR: unexpected value for norm_method")
} else {
    logger_ts("Note: not focusing solely on metabolic genes")
    pseudo_norm <- cpm(dgelist_metab, log = FALSE)
}

logger_ts("Outputting Expression Matrix")
write.table(
    pseudo_norm,
    file = output_path("pseudo_matrix"),
    sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE
)

logger_ts("Outputting Metadata")
write.table(
    pseudo_meta,
    output_path("pseudo_metadata"),
    sep = "\t", col.names = TRUE, row.names = TRUE, quote = TRUE
)
output_json(summarize_discrete_metadata(pseudo_meta), 'discrete_metadata_summary')

logger_ts("Outputting cells and samples lists")
output_json(celltypes, 'cell_types')
output_json(celltypes, 'samples')

logger_ts("DONE")
