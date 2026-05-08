
# ============================================================
# IGF2BP2 PROJECT – Rhythmicity analysis
# ============================================================

# Author: Shriyansh Srivastava
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script performs Gene ontology analysis using the Clusterprofiler package

# Input:
# Gene expression matrix from publicly available datasets
#
# Output:
# Rhythmicity of genes assigned by the Metacycle algorithms



# ============================================================
# FIGURE S3 - The outputs from this was used for plotting Figure S3
# along with the results from RAIN
# ============================================================




setwd("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/rhythmicity_analyses/metacycle/")
getwd()


library(MetaCycle)
library(dplyr)
library(ggplot2)


### datasets used in the current study for rhythmicity analysis ###
 
# # GSE197726 - doi:10.1016/j.celrep.2023.112588####
# # "WT
# ZT00 – 3
# ZT04 – 3
# ZT08 – 3
# ZT12 – 3
# ZT16 – 3
# ZT20 – 3
# Total WT – 18
# 
# KO
# ZT00 – 3
# ZT04 – 3
# ZT08 – 3
# ZT12 – 3
# ZT16 – 3
# ZT20 – 3
# Total KO – 18
# 
# mRE
# ZT00 – 3
# ZT04 – 3
# ZT08 – 3
# ZT12 – 3
# ZT16 – 3
# ZT20 – 3
# Total mRE – 18
# 
# Total - 54"

# ---------- read data ----------

GSE197726_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/output_text_files/07_filtered_counts/gse197726_normalized_counts.tsv",
  stringsAsFactors = FALSE,
  row.names = 1
)

# ---------- experiment design ----------

timepoints <- rep(c(0,4,8,12,16,20), each = 3)

conditions <- c("WT","KO","mRE")


# ---------- loop for metacycle run ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(paste0("^", cond), colnames(GSE197726_countData), value = TRUE)
  mat  <- GSE197726_countData[, cols]
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  # create directories
  infile_dir  <- file.path("dataset/GSE197726/metacycle_infile")
  output_dir  <- file.path("dataset/GSE197726/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(infile_dir,
                      paste0("GSE197726_", cond, "_metacycle_infile.txt"))
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
  
}





# # GSE194106 - doi:10.1038/s42255-023-00826-7####
# "DRF
# ZT00 – 4
# ZT04 – 4
# ZT08 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# ZT24 – 4
# Total DRF – 28
# 
# NRF
# ZT00 – 4
# ZT04 – 4
# ZT08 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# ZT24 – 4
# Total NRF – 28
# 
# Total - 56"


# ---------- read data ----------

GSE194106_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Mouse/GSE194106/GSE194106_RNAseq-muscle-tpm.txt",
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

# igf2bp2 319765

# ---------- experiment design ----------

conditions <- c("GasD","GasN")

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(cond, colnames(GSE194106_countData), value = TRUE)
  mat <- GSE194106_countData[, cols]
  
  # extract time (GasD00 → 0)
  times <- as.numeric(sub(".*Gas[DN](\\d+).*", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  # ordered timepoints
  timepoints <- as.numeric(sub(".*Gas[DN](\\d+).*", "\\1", colnames(mat)))
  
  # QC checks
  colnames(mat)[1:20]
  table(timepoints)
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  # directories
  infile_dir <- file.path("dataset/GSE194106/metacycle_infile")
  output_dir <- file.path("dataset/GSE194106/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE194106_", cond, "_metacycle_infile.txt")
  )
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
}



# # GSE273878 - doi:10.1016/j.celrep.2025.115689####
# "WT_CTL (wild-type sham/control)
# 00:00 — 2
# 04:00 — 2
# 08:00 — 2
# 12:00 — 2
# 16:00 — 2
# 20:00 — 2
# Total WT_CTL — 12
# 
# WT_KPC
# 00:00 — 2
# 04:00 — 2
# 08:00 — 2
# 12:00 — 2
# 16:00 — 2
# 20:00 — 2
# Total WT_KPC — 12
# 
# KO_CTL (FoxP1 skeletal muscle KO sham/control)
# 00:00 — 2
# 04:00 — 2
# 08:00 — 2
# 12:00 — 2
# 16:00 — 2
# 20:00 — 2
# Total KO_CTL — 12
# 
# KO_KPC (FoxP1 skeletal muscle KO + KPC)
# 00:00 — 2
# 04:00 — 2
# 08:00 — 2
# 12:00 — 2
# 16:00 — 2
# 20:00 — 2
# Total KO_KPC — 12
# 
# Total - 48"


# ---------- read data ----------

GSE273878_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Mouse/GSE273878/GSE273878_gene_fpkm.xls",
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

# remove annotation columns if present
GSE273878_countData <- GSE273878_countData[,-(49:57)]

# ---------- experiment design ----------

conditions <- c("SW","KW","SM","KM") #SW = WT_CTL, KW = WT_KPC, SM = KO_CTL, KM = KO_KPC

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(paste0("^",cond), colnames(GSE273878_countData), value = TRUE)
  mat <- GSE273878_countData[, cols]
  
  # extract time (SW04_1 → 4)
  times <- as.numeric(sub(".*([0-9]{2})_.*", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  # ordered timepoints
  timepoints <- as.numeric(sub(".*([0-9]{2})_.*", "\\1", colnames(mat)))
  
  # QC checks
  colnames(mat)
  table(timepoints)
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  # directories
  infile_dir <- file.path("dataset/GSE273878/metacycle_infile")
  output_dir <- file.path("dataset/GSE273878/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE273878_", cond, "_metacycle_infile.txt")
  )
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
}




# # GSE308276 - doi:10.1177/07487304251386926####
# "WT
# ZT0 – 3
# ZT4 – 3
# ZT8 – 3
# ZT12 – 3
# ZT16 – 3
# ZT20 – 3
# Total WT – 18
# 
# Bmal1hep-/-
# ZT0 – 3
# ZT4 – 3
# ZT8 – 3
# ZT12 – 3
# ZT16 – 3
# ZT20 – 3
# Total Bmal1hep-/- – 18
# 
# Total - 36"

# ---------- read data ----------

library(readxl)
GSE308276_countData <- read_xlsx(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Mouse/GSE308276/GSE308276_FPKM_mus_from_WT_and_Bmal1LKO_6tp.xlsx")

colnames(GSE308276_countData)


GSE308276_countData <- as.data.frame(GSE308276_countData)

# clean column names (remove Excel prefixes like "B1: ")
colnames(GSE308276_countData) <- gsub("B[0-9]+: ", "", colnames(GSE308276_countData))

# set gene IDs
rownames(GSE308276_countData) <- GSE308276_countData$GeneSym
GSE308276_countData <- GSE308276_countData[,-1]

# ---------- experiment design ----------

conditions <- c("WT","KO")

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(cond, colnames(GSE308276_countData), value = TRUE)
  mat <- GSE308276_countData[, cols]
  
  # extract time (ZT format)
  times <- as.numeric(sub(".*ZT([0-9]+).*", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  # ordered timepoints
  timepoints <- as.numeric(sub(".*ZT([0-9]+).*", "\\1", colnames(mat)))
  
  # QC checks
  colnames(mat)
  table(timepoints)
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  # directories
  infile_dir <- file.path("dataset/GSE308276/metacycle_infile")
  output_dir <- file.path("dataset/GSE308276/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE308276_", cond, "_metacycle_infile.txt")
  )
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
}




# # GSE107787 - doi:10.1016/j.celrep.2018.11.077####
# "Liver_FAST
# 00:00 — 3
# 04:00 — 3
# 08:00 — 3
# 12:00 — 3
# 16:00 — 3
# 20:00 — 3
# Total Liver_FAST — 18
# 
# Liver_FED
# 00:00 — 3
# 04:00 — 3
# 08:00 — 3
# 12:00 — 3
# 16:00 — 3
# 20:00 — 3
# Total Liver_FED — 18
# 
# Muscle_FAST
# 00:00 — 3
# 04:00 — 3
# 08:00 — 3
# 12:00 — 3
# 16:00 — 3
# 20:00 — 3
# Total Muscle_FAST — 18
# 
# Muscle_FED
# 00:00 — 3
# 04:00 — 3
# 08:00 — 3
# 12:00 — 3
# 16:00 — 3
# 20:00 — 3
# Total Muscle_FED — 18
# 
# Total - 72"


# ---------- read data ----------

GSE107787_countData <- read_xlsx(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Mouse/GSE107787/GSE107787_RNAseq_expression_results.xlsx")

GSE107787_countData <- as.data.frame(GSE107787_countData)

# move gene names to rownames
rownames(GSE107787_countData) <- GSE107787_countData$Gene
GSE107787_countData <- GSE107787_countData[,-1]

# ---------- experiment design ----------

conditions <- c("Liver_FAST","Liver_FED","Muscle_FAST","Muscle_FED")

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(cond, colnames(GSE107787_countData), value = TRUE)
  mat <- GSE107787_countData[, cols]
  
  # extract time
  times <- as.numeric(sub(".*ZT([0-9]+).*", "\\1", cols))
  
  # reorder columns
  mat <- mat[, order(times)]
  
  # ordered timepoints
  timepoints <- as.numeric(sub(".*ZT([0-9]+).*", "\\1", colnames(mat)))
  
  # QC
  colnames(mat)
  table(timepoints)
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  infile_dir <- file.path("dataset/GSE107787/metacycle_infile")
  output_dir <- file.path("dataset/GSE107787/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE107787_", cond, "_metacycle_infile.txt")
  )
  
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
}




# # GSE195724 - doi:10.1126/science.adj8533 ###
# "ZT0 – 3
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – WT_10w_ALF – 23
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 3
# ZT12 – 4
# ZT16 – 4
# ZT20 – 3
# 
# Total – WT_26w_ALF – 22
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 5
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – Mu-RE_10w_ALF – 25
# 
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 3
# ZT20 – 5
# 
# Total – Mu-RE_26w_ALF – 24
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 3
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – Br-RE_10w_ALF – 23
# 
# ZT0 – 4
# ZT4 – 3
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 3
# 
# Total – Br-RE_26w_ALF – 22
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – RE-RE_10w_ALF – 24
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – RE-RE_26w_ALF – 24
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 3
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – KO_10w_ALF – 23
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 5
# ZT20 – 3
# 
# Total – KO_26w_ALF – 24
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – WT_26w_TRF – 24
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 3
# ZT16 – 4
# ZT20 – 4
# 
# Total – Mu-RE_26w_TRF – 23
# 
# ZT0 – 4
# ZT4 – 3
# ZT8 – 4
# ZT12 – 5
# ZT16 – 3
# ZT20 – 4
# 
# Total – KO_26w_TRF – 23
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – Old_96w_ALF – 24
# 
# ZT0 – 4
# ZT4 – 4
# ZT8 – 4
# ZT12 – 4
# ZT16 – 4
# ZT20 – 4
# 
# Total – Old_96w_TRF – 24
# 
# Total - 352"


# ---------- read data ----------

GSE195724_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Mouse/GSE195724_RNAseq_gene_expression_log2FPKM.txt",
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)


# ---------- experiment design ----------

conditions <- c(
  "WT_10w_ALF","WT_26w_ALF",
  "Mu-RE_10w_ALF","Mu-RE_26w_ALF",
  "Br-RE_10w_ALF","Br-RE_26w_ALF",
  "RE-RE_10w_ALF","RE-RE_26w_ALF",
  "KO_10w_ALF","KO_26w_ALF",
  "WT_26w_TRF","Mu-RE_26w_TRF",
  "KO_26w_TRF",
  "Old_96w_ALF","Old_96w_TRF"
)


# ---------- loop for MetaCycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  cols <- grep(cond, colnames(GSE195724_countData), value = TRUE)
  
  mat <- GSE195724_countData[, cols]
  
  # extract time from column names
  times <- as.numeric(sub(".*ZT([0-9]+).*", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  # recompute timepoints AFTER sorting
  timepoints <- as.numeric(sub(".*ZT([0-9]+).*", "\\1", colnames(mat)))
  print(table(timepoints))
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  infile_dir <- file.path("dataset/GSE195724/metacycle_infile")
  output_dir <- file.path("dataset/GSE195724/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE195724_", cond, "_metacycle_infile.txt")
  )
  
  write.table(
    meta_df,
    infile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30
  )
  
  gc()
  
}




# # GSE182117 - doi:10.1126/sciadv.abi9654 ####
# "NGT_CTL
# 12:00 — 7
# 18:00 — 6
# 24:00 — 6
# 30:00 — 7
# 36:00 — 7
# 42:00 — 6
# 48:00 — 7
# 54:00 — 7
# Total NGT_CTL - 53
# 
# NGT_HGI
# 12:00 — 7
# 18:00 — 6
# 24:00 — 6
# 30:00 — 7
# 36:00 — 7
# 42:00 — 7
# 48:00 — 7
# 54:00 — 6
# Total NGT_HGI - 53
# 
# T2D_CTL
# 12:00 — 5
# 18:00 — 5
# 24:00 — 5
# 30:00 — 5
# 36:00 — 5
# 42:00 — 5
# 48:00 — 5
# 54:00 — 5
# Total T2D_CTL - 40
# 
# T2D_HGI
# 12:00 — 5
# 18:00 — 5
# 24:00 — 5
# 30:00 — 5
# 36:00 — 5
# 42:00 — 5
# 48:00 — 5
# 54:00 — 5
# Total T2D_HGI - 40
# 
# Total - 186"
# 

# ---------- read data ----------

GSE182117_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Human/GSE182117/GSE182117_logCPM_batchCorrected.tsv",
  stringsAsFactors = FALSE,
  row.names = 1
)

# ---------- experiment design ----------

conditions <- c("NGT_CTL","NGT_HGI","T2D_CTL","T2D_HGI")

timepoints_list <- list(
  
  NGT_CTL = c(rep(12,7), rep(18,6), rep(24,6), rep(30,7),
              rep(36,7), rep(42,6), rep(48,7), rep(54,7)),
  
  NGT_HGI = c(rep(12,7), rep(18,6), rep(24,6), rep(30,7),
              rep(36,7), rep(42,7), rep(48,7), rep(54,6)),
  
  T2D_CTL = c(rep(12,5), rep(18,5), rep(24,5), rep(30,5),
              rep(36,5), rep(42,5), rep(48,5), rep(54,5)),
  
  T2D_HGI = c(rep(12,5), rep(18,5), rep(24,5), rep(30,5),
              rep(36,5), rep(42,5), rep(48,5), rep(54,5))
)

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(paste0("^", cond), colnames(GSE182117_countData), value = TRUE)
  mat <- GSE182117_countData[, cols]
  
  # extract time from column names
  times <- as.numeric(sub(".*_(\\d+)h_.*", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  timepoints <- as.numeric(sub(".*_(\\d+)h_.*", "\\1", colnames(mat)))
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  colnames(mat)[1:53]
  table(times)
  
  # directories
  infile_dir <- file.path("dataset/GSE182117/metacycle_infile")
  output_dir <- file.path("dataset/GSE182117/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(infile_dir,
                      paste0("GSE182117_", cond, "_metacycle_infile.txt"))
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
  
}



# # GSE109825 - doi:10.7554/eLife.34114 ####
# "siControl
# 00:00 — 2
# 02:00 — 2
# 04:00 — 2
# 06:00 — 2
# 08:00 — 2
# 10:00 — 2
# 12:00 — 2
# 14:00 — 2
# 16:00 — 2
# 18:00 — 2
# 20:00 — 2
# 22:00 — 2
# 24:00 — 2
# 26:00 — 2
# 28:00 — 2
# 30:00 — 2
# 32:00 — 2
# 34:00 — 2
# 36:00 — 2
# 38:00 — 2
# 40:00 — 2
# 42:00 — 2
# 44:00 — 2
# 46:00 — 2
# 48:00 — 2
# Total siControl — 50
# 
# siClock
# 00:00 — 2
# 02:00 — 2
# 04:00 — 2
# 06:00 — 2
# 08:00 — 2
# 10:00 — 2
# 12:00 — 2
# 14:00 — 2
# 16:00 — 2
# 18:00 — 2
# 20:00 — 2
# 22:00 — 2
# 24:00 — 2
# 26:00 — 2
# 28:00 — 2
# 30:00 — 2
# 32:00 — 2
# 34:00 — 2
# 36:00 — 2
# 38:00 — 2
# 40:00 — 2
# 42:00 — 2
# 44:00 — 2
# 46:00 — 2
# 48:00 — 2
# Total siClock — 50
# 
# Total - 100"


# ---------- read data ----------

GSE109825_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Human/GSE109825/GSE109825_renamed.tsv", 
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names=1
)
#igf2bp2 ENSG00000073792.11
GSE109825_countData <- GSE109825_countData[,-(1:2)]

dput(colnames(GSE109825_countData)[1:10])
unique(gsub(".*_(siControl|siClock).*", "\\1", colnames(GSE109825_countData)))

# ---------- experiment design ----------

conditions <- c("siControl","siClock")

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select columns
  cols <- grep(cond, colnames(GSE109825_countData), value = TRUE)
  mat <- GSE109825_countData[, cols]
  mat
  
  # extract time from column names
  times <- as.numeric(sub(".*ZT([0-9]+)", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  # extract ordered timepoints
  timepoints <- as.numeric(sub(".*ZT([0-9]+)", "\\1", colnames(mat)))
  
  colnames(mat)[1:50]
  table(times)
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  # directories
  infile_dir <- file.path("dataset/GSE109825/metacycle_infile")
  output_dir <- file.path("dataset/GSE109825/metacycle_output", cond)
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE109825_", cond, "_metacycle_infile.txt")
  )
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
}






# # GSE108539 - doi:10.7554/eLife.34114####
# # "00:00 — 9
# 04:00 — 9
# 08:00 — 9
# 12:00 — 10
# 16:00 — 10
# 20:00 — 10
# 
# Total - 57"

# ---------- read data ----------

GSE108539_countData <- read.delim(
  "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/circadian_datasets_for_rhythmicity_analysis/Human/GSE108539/GSE108539_norm_reads_count.txt",
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

# igf2bp2 ENSG00000073792

# remove annotation columns if present
GSE108539_countData <- GSE108539_countData[,-(1:57)]

# ---------- experiment design ----------

# single condition dataset
conditions <- c("all_samples")

# ---------- loop for metacycle ----------

for (cond in conditions) {
  
  message("Processing ", cond)
  
  # select all columns
  cols <- colnames(GSE108539_countData)
  mat <- GSE108539_countData[, cols]
  
  # extract time from column names (example: sample_04h)
  times <- as.numeric(sub(".*_(\\d+)\\.\\d+", "\\1", cols))
  
  # reorder columns by time
  mat <- mat[, order(times)]
  
  # extract ordered timepoints
  timepoints <- as.numeric(sub(".*_(\\d+)\\.\\d+", "\\1", colnames(mat)))
  
  # quick QC checks
  colnames(mat)[1:20]
  table(timepoints)
  
  # prepare MetaCycle input
  meta_df <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
  
  # directories
  infile_dir <- file.path("dataset/GSE108539/metacycle_infile")
  output_dir <- file.path("dataset/GSE108539/metacycle_output")
  
  dir.create(infile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  infile <- file.path(
    infile_dir,
    paste0("GSE108539_metacycle_infile.txt")
  )
  
  # write infile
  write.table(meta_df,
              infile,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  # run MetaCycle
  meta2d(
    infile = infile,
    outdir = output_dir,
    filestyle = "txt",
    timepoints = timepoints,
    minper = 20,
    maxper = 28,
    cycMethod = c("ARS","JTK","LS"),
    analysisStrategy = "auto",
    outputFile = TRUE,
    outIntegration = "both",
    adjustPhase = "predictedPer",
    combinePvalue = "fisher",
    weightedPerPha = FALSE,
    ARSmle = "auto",
    ARSdefaultPer = 24,
    outRawData = FALSE,
    releaseNote = TRUE,
    outSymbol = "output",
    parallelize = TRUE,
    nCores = 30,
    inDF = NULL
  )
  
  gc()
}




save.image(file = "script/metacycle_analysis.RData")
