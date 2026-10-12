#!/usr/bin/env Rscript
# Per-sample QC: filter low-quality cells, then remove doublets with scDblFinder.
# Input: a 10x Cell Ranger "filtered_feature_bc_matrix" directory.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scDblFinder)
})

option_list <- list(
  make_option("--sample",       type = "character"),
  make_option("--condition",    type = "character"),
  make_option("--donor",        type = "character"),
  make_option("--matrix",       type = "character"),
  make_option("--min_features", type = "integer",   default = 200),
  make_option("--max_mito",     type = "double",    default = 10),
  make_option("--mito_pattern", type = "character", default = "^MT-")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$sample), !is.null(opt$matrix))

set.seed(1)

counts <- Read10X(data.dir = opt$matrix)
if (is.list(counts)) counts <- counts[["Gene Expression"]]

obj <- CreateSeuratObject(counts = counts, project = opt$sample, min.cells = 3)
obj$sample    <- opt$sample
obj$condition <- opt$condition
obj$donor     <- opt$donor
# Sanity check: a mito pattern that matches nothing usually means the wrong species
# or annotation (e.g. mouse gene names such as 'mt-Nd1' with a human '^MT-' pattern).
n_mito_genes <- sum(grepl(opt$mito_pattern, rownames(obj)))
if (n_mito_genes == 0) {
  stop("No genes match mito_pattern '", opt$mito_pattern, "' in sample ", opt$sample,
       ". Check the species/annotation in features.tsv (first genes: ",
       paste(head(rownames(obj), 5), collapse = ", "), ").")
}
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = opt$mito_pattern)
n_input <- ncol(obj)

# 1. Basic QC filters
keep <- obj$nFeature_RNA >= opt$min_features & obj$percent.mt <= opt$max_mito
obj <- subset(obj, cells = colnames(obj)[keep])
n_basic <- ncol(obj)

# 2. Doublet detection on the QC-filtered cells
sce <- scDblFinder(GetAssayData(obj, layer = "counts"))
obj$doublet_class <- as.character(sce$scDblFinder.class)
n_doublets <- sum(obj$doublet_class == "doublet")
obj <- subset(obj, cells = colnames(obj)[obj$doublet_class == "singlet"])

metrics <- data.frame(
  sample              = opt$sample,
  condition           = opt$condition,
  donor               = opt$donor,
  cells_input         = n_input,
  cells_after_qc      = n_basic,
  doublets_removed    = n_doublets,
  cells_final         = ncol(obj),
  median_genes        = median(obj$nFeature_RNA),
  median_umis         = median(obj$nCount_RNA),
  median_percent_mito = round(median(obj$percent.mt), 2)
)

write.table(metrics, paste0(opt$sample, ".qc_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(obj, paste0(opt$sample, ".qc.rds"))
writeLines(capture.output(sessionInfo()), paste0(opt$sample, ".sessionInfo.txt"))
