############################################################
# IGF2BP2 PROJECT – Differential Expression Analysis (DESeq2)
############################################################
# Author: Vishnu A M
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
# This script performs differential gene expression analysis
# between WT and KO samples using DESeq2.
#
# Input:
# - Raw count matrix (dds75_filtered_raw_counts.tsv). The input file contain 13918 genes
# that was used in the WGCNA analysis. The corresponding dataset's GEO identifier: GSE197726
#
# Output:
# - Differential expression results (CSV)
#
# Notes:
# - Timepoint is included as a covariate
# - Volcano plotting handled separately (Python)

##############################
# Paths
##############################
input_counts_file <- "path/to/input/dds75_filtered_raw_counts.tsv"
output_dir <- "path/to/output"
output_file <- file.path(output_dir, "gse197726_salmon_deseq2_results.csv")

##############################
# Libraries
##############################
library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tibble)
library(EnhancedVolcano)

##############################
# Load data
##############################
counts <- read.table(
  input_counts_file,
  header = TRUE,
  row.names = 1
)

print(colnames(counts))

##############################
# Subset WT and KO
##############################
counts_wt_ko <- counts[, grepl("^(WT|KO)\\.", colnames(counts))]
counts_wt_ko <- round(counts_wt_ko)

##############################
# Metadata
##############################
sample_names <- colnames(counts_wt_ko)

condition <- factor(
  ifelse(startsWith(sample_names, "WT."), "WT", "KO"),
  levels = c("WT", "KO")
)

timepoint <- factor(
  sub(".*ZT\\.([0-9]+).*", "ZT_\\1", sample_names)
)

coldata <- data.frame(
  row.names = sample_names,
  condition = condition,
  timepoint = timepoint
)

##############################
# DESeq2
##############################
dds <- DESeqDataSetFromMatrix(
  countData = counts_wt_ko,
  colData = coldata,
  design = ~ timepoint + condition
)

dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)

res <- results(dds)
summary(res)

##############################
# Output directory
##############################
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

##############################
# Process results
##############################
res_df <- as.data.frame(res)

rownames(res_df) <- toupper(rownames(res_df))

res_df <- rownames_to_column(res_df, var = "gene")

##############################
# Annotation
##############################
res_df$symbol <- mapIds(
  org.Mm.eg.db,
  keys = res_df$gene,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

##############################
# Save
##############################
write.csv(res_df, output_file, row.names = FALSE)

############################################################
# Optional plotting (disabled)
############################################################
# genes_of_interest <- c("Igf2bp2", "Nr1d2", "Nr1d1")

# plot <- EnhancedVolcano(
#   res_df,
#   x = "log2FoldChange",
#   y = "padj",
#   lab = res_df$symbol,
#   pCutoff = 0.01,
#   FCcutoff = 1,
#   selectLab = genes_of_interest,
#   drawConnectors = TRUE,
#   boxedLabels = TRUE,
#   max.overlaps = Inf,
#   title = "Wildtype vs Knockout",
#   xlim = c(-5, 5)
# )

# print(plot)