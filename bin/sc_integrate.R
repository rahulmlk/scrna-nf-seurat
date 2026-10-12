#!/usr/bin/env Rscript
# Merge QC'd samples, normalise, run PCA, correct batch with Harmony, cluster, embed.
# Usage: sc_integrate.R [options] sample1.qc.rds sample2.qc.rds ...

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

option_list <- list(
  make_option("--batch_var",  type = "character", default = "donor"),
  make_option("--n_pcs",      type = "integer",   default = 30),
  make_option("--resolution", type = "double",    default = 0.5)
)
parsed <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
opt   <- parsed$options
files <- parsed$args
stopifnot(length(files) >= 2)

set.seed(1)

objs <- lapply(files, readRDS)
ids  <- vapply(objs, function(o) as.character(o$sample[1]), character(1))

merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = ids)
merged <- JoinLayers(merged)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = opt$n_pcs, verbose = FALSE)

# Correct only for a technical/batch variable. Never correct for the condition you
# want to test, or Harmony can remove the biology. If every level of the batch
# variable belongs to a single condition (e.g. one sample per donor, donors split
# between conditions), correcting for it is the same as correcting for condition.
use_harmony <- !(tolower(opt$batch_var) %in% c("none", "null", "na", ""))
reduction   <- "pca"

if (use_harmony) {
  md <- merged@meta.data
  stopifnot(opt$batch_var %in% colnames(md))
  tab    <- table(md[[opt$batch_var]], md$condition)
  nested <- all(rowSums(tab > 0) == 1) && ncol(tab) > 1
  if (nested) {
    stop("batch_var '", opt$batch_var, "' is nested within condition, so Harmony would ",
         "remove the condition effect. Use --batch_var none, or a technical batch ",
         "variable that spans conditions (e.g. sequencing run, chemistry).")
  }
  merged    <- RunHarmony(merged, group.by.vars = opt$batch_var, reduction.use = "pca")
  reduction <- "harmony"
}

dims <- seq_len(opt$n_pcs)
merged <- FindNeighbors(merged, reduction = reduction, dims = dims)
merged <- FindClusters(merged, resolution = opt$resolution)
merged <- RunUMAP(merged, reduction = reduction, dims = dims)

saveRDS(merged, "integrated.rds")

counts_tbl <- as.data.frame(table(cluster = merged$seurat_clusters, sample = merged$sample))
write.table(counts_tbl, "cluster_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

plot_vars <- c("seurat_clusters", "condition", "sample")
if (use_harmony) plot_vars <- unique(c(plot_vars, opt$batch_var))
p <- DimPlot(merged, reduction = "umap", group.by = plot_vars, combine = TRUE)
ggsave("umap.png", p, width = 5 * length(plot_vars), height = 5, dpi = 150)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
