############################################################
# IGF2BP2 PROJECT – tximport (Salmon → Gene Counts)
############################################################

# Author: Vishnu A M
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script gets the gene level counts from Salmon quant files

# Input:
# Salmon quant files
#
# Output:
# Gene level counts file as a tsv
#


##############################
# Paths
##############################
quants_dir <- "path/to/salmon_quant"
output_file <- "path/to/output/salmon_gene_counts.tsv"

##############################
# Libraries
##############################
library(tximport)
library(ensembldb)
library(AnnotationHub)
library(GEOquery)

##############################
# 1. Transcript → Gene mapping
##############################
hub <- AnnotationHub()

ensdb_query <- query(hub, c("EnsDb", "Mus musculus", "113"))
ensdb_113 <- ensdb_query[["AH119358"]]

tx_data <- transcripts(ensdb_113, return.type = "DataFrame")
tx2gene <- tx_data[, c("tx_id", "gene_id")]

##############################
# 2. Load Salmon quant files
##############################
quant_files <- list.files(
  quants_dir,
  pattern = "quant.sf$",
  recursive = TRUE,
  full.names = TRUE
)

quant_dirs <- list.files(
  quants_dir,
  pattern = "_quant$",
  full.names = TRUE
)

sample_names <- gsub("_quant$", "", basename(quant_dirs))
names(quant_files) <- sample_names

##############################
# 3. tximport
##############################
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

raw_counts <- txi$counts
raw_counts_df <- data.frame(raw_counts)

##############################
# 4. Metadata (GEO)
##############################
geo_id <- "GSE197726"

gse_metadata <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse_metadata[[1]]))
phenoData <- phenoData[c(1, 2, 11, 12)]

##############################
# 5. Assign sample names
##############################
colnames(raw_counts_df) <- phenoData$title

##############################
# 6. Save output
##############################
write.table(
  raw_counts_df,
  file = output_file,
  sep = "\t",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)