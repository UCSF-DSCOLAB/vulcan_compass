suppressPackageStartupMessages({
    library(Seurat)
})

### Parameters
scd_rds <- input_path('scd_rds')
min_cells <- params['min_cells']
samp_col <- params['sample_id_column']
ct_col <- params['cell_type_column']
norm_method <- params['norm_method']

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

# Metadata, per-cell
meta <- sobj@meta.data

# Counts - delog
norm <- GetAssayData(sobj, assay = 'RNA', slot = 'data')
rm(sobj)
gc(verbose=FALSE)
delog_norm <- norm
delog_norm@x <- exp(delog_norm@x)
stopifnot(identical(delog_norm[5,6], exp(norm[5,6])))
rm(norm)
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
### Function for summarizing metadat for later UIs
summarize_discrete_metadata <- function(metadata) {
    
    output <- list()
    for (i in getMetas(scdata)) {
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
# Donors/Samples
cell_samples <- meta[,samp_col]
samples <- unique(cell_samples)
# Tracking skippage due to min_cells
skipped_sample <- NULL
skipped_cell <- NULL
skipped_i <- NULL
i <- 0
# Tracking metadata skippage
meta_ignored <- NULL

if (norm_method=="early__depth_norm_to_metabolic_genes") {
    logger_ts("Normalizing counts per cell toward metabolic genes")
    logger_ts("NOT YET IMPLEMENTED")
}

logger_ts("Starting Pseudobulking")
pseudo_mat <- NULL
pseudo_meta <- NULL
for (ct in celltypes) {
    logger_ts("\tworking on ", ct)
    is_cell <- cell_annots==ct
    for (samp in samples) {
        i <- i + 1
        is_set <- is_cell & cell_samples==samp
        if (sum(is_set) > min_cells) {
            # Trim matrix and make dense
            delog_norm_set <- as.matrix(delog_norm[,is_set])
            meta_set <- meta[is_set,]
            this_pseudo_mat <- apply(delog_norm_set, 1, mean)
            this_pseudo_meta <- sapply(1:ncol(meta_set), function(meta_i) {
                meta_collapse(meta_set[,meta_i], names(meta_set)[meta_i])
            })
            this_name <- paste0(samp, "__", ct)
            # Add column
            if (identical(pseudo_mat, NULL)) {
                pseudo_mat <- data.frame(this_pseudo_mat)
                colnames(pseudo_mat) <- this_name
                pseudo_meta <- meta[1,]
                pseudo_meta[1,] <- this_pseudo_meta
                rownames(pseudo_meta) <- this_name
            } else {
                pseudo_mat <- cbind(pseudo_mat, this_pseudo_mat)
                colnames(pseudo_mat)[ncol(pseudo_mat)] <- this_name
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
rownames(pseudo_mat) <- rownames(delog_norm)

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

logger_ts("Outputting Expression Matrix")
write.table(
    pseudo_mat,
    file = output_path(pseudo_matrix),
    sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE
)

logger_ts("Outputting Metadata")
write.table(
    pseudo_meta,
    output_path(pseudo_metadata),
    sep = "\t", col.names = TRUE, row.names = TRUE, quote = TRUE
)
output_json(summarize_discrete_metadata(pseudo_meta), 'discrete_metadata_summary')

logger_ts("Outputting cells and samples lists")
output_json(celltypes, 'cell_types')
output_json(celltypes, 'samples')

logger_ts("DONE")
