############################################################
# IGF2BP2 PROJECT – WGCNA
############################################################

# Author: Vishnu A M
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script runs WGCNA on the Tximport generated raw gene expression matrix

# Input:
# Raw gene expression matrix
#
# ### Output Figures ###
#
# Figure 1:
# Figure 1A - Hierarchical clustering
# Figure 1B - Principal Component Analysis
# Figure 1C - No. of. Genes per Module
# Figure 1D - Module-trait correlation 
#
# Figure S1:
# Figure S1
#
#
# ### Output Text files ###
# Supplementary file 1 - Gene module list
# Supplementary file 3 - GS and MM
# 
#


# ============================
# Libraries
# ============================
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(dendextend)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)

# Enable multithreading
allowWGCNAThreads()

#set the working directory
setwd("/set/path")

#save the workspace image
save.image("/set path/.RData")


# 1. Load the data
#load the workspace image
load("/set path/.RData")


# 1. Load Data and prepare the data ------------------------------------------------Differential_gene_expression.R

data <- read.delim("/path/to/gene/expression/matrix", header = TRUE)

print(colnames(data))


#get metadata
geo_id <- "GSE197726"
gse_metadata <- getGEO (geo_id, GSEMatrix = TRUE)
phenoData<- pData (phenoData(gse_metadata[[1]]))
head (phenoData)
phenoData<- phenoData[c(1,2,11,12)]

# change the inconsistencies in PhenoData title column
phenoData$title <- gsub("_", "-", phenoData$title)
phenoData$title <- gsub(" ", "-", phenoData$title)


# Replace column names of data (excluding first column)
colnames(data)[-1] <- phenoData$title
colnames(data)

# prepare data
data[1:10,1:10]
print(colnames(data))

data <- data %>%
  pivot_longer(cols = -"ensembl_id", names_to = "samples", values_to = "counts") %>%
  inner_join(phenoData, by = c("samples" = "title")) %>%
  dplyr::select(1,2,3) %>%
  pivot_wider(names_from = samples, values_from = counts, names_sort = FALSE) %>%
  column_to_rownames(var = "ensembl_id")

# 2. QC - outlier detection ------------------------------------------------

# detect outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData


# fixing column names in colData
colnames(colData) <- c("title", "geo_accession", "genotype", "timepoint")

#cleaning up colData
colData <- colData %>%
  mutate(
    genotype = gsub("genotype: ", "", genotype),
    timepoint = gsub("time point: ", "", timepoint)
  )


#change the rownames
rownames(colData) <- colData$title

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data))
all(rownames(colData) == colnames(data))

#round the values of the data variable as deseq expects integers
countData<-round(data)


# create dds
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ 1) # not specifying model


## remove all genes with counts < 15 in more than 75% of samples (54*0.75=40.5 )
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 41,]
nrow(dds75) # 13918 genes


counts_df <- counts(dds75)

## This output file is used for running Differential gene expression. The code corresponding to that
# is named as "Differential_gene_expression.R"
write.table(
  counts_df,
  file = "set/path/dds75_filtered_raw_counts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

norm.counts_t<- t(norm.counts)
# save norm.counts
write.table(
  norm.counts_t,
  file = "set/path/output",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
# =======================================================
# FIGURE 1A
# =======================================================


# Perform Hierarchichal clustering and PCA

###### detect outlier samples - hierarchical clustering - method 1

# Define output path
outpath_hclust <- "set/output/path"


htree <- hclust(dist(norm.counts), method = "average")
dend <- as.dendrogram(htree)

# Extract groups and map colors
labels_list <- labels(dend)
group <- gsub("-.*", "", labels_list)
group_colors <- c("WT" = "orange", "KO" = "darkblue", "mRE" = "brown")

# --- AESTHETICS ---
labels_colors(dend) <- group_colors[group]
labels_cex(dend)    <- 0.5  # Slightly smaller to accommodate bold width


# 3. Save the Plot (Set to 5x3 as requested)
tiff(
  filename = file.path(outpath_hclust, "01_htree_final_bold.tiff"),
  width = 6,
  height = 3,
  units = "in",
  res = 300,
  bg = "transparent"
)

# 4. Plotting Parameters
par(
  mar = c(1, 4, 1, 0), # Increased bottom margin for vertical bold labels
  family = "sans",
  las = 2              # Rotates labels to be vertical (essential for 5x3)
)

# 5. Execute Plot
plot(
  hang.dendrogram(dend, hang = 0.1), 
  main = NULL, 
  ylab = "Height"
)

# 6. Legend
legend(
  "topright",
  legend = names(group_colors),
  fill = group_colors,
  border = "black",
  cex = 0.6,           # Smaller legend to fit 5x3 scale
  bty = "n",
  inset = c(0.01, 0.01)
)

dev.off()


# =======================================================
# FIGURE 1B
# =======================================================

###### detect outlier samples - Principal Component Analysis (PCA) - method 2

## -------------------------------
## 1. PCA computation
## -------------------------------
pca <- prcomp(norm.counts)

# % variance explained
pca.var.percent <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

# PCA scores
pca.dat <- as.data.frame(pca$x)
pca.dat$samples <- rownames(pca.dat)


## -------------------------------
## 2. Automatically assign genotype 
##    based on sample names
## -------------------------------
# Detect WT, KO, mRE directly from sample ID
pca.dat$genotype <- case_when(
  grepl("^WT",  pca.dat$samples, ignore.case = TRUE)  ~ "WT",
  grepl("^KO",  pca.dat$samples, ignore.case = TRUE)  ~ "KO",
  grepl("^mRE", pca.dat$samples, ignore.case = TRUE)  ~ "mRE",
  TRUE                                                ~ "Unknown"
)

# Check results
print(table(pca.dat$genotype))

## Save PCA scores
write.csv(pca.dat,
          file = "../01.1_pca/PCA_scores.csv",
          row.names = FALSE)

## Save variance explained
write.csv(data.frame(PC = paste0("PC", 1:length(pca.var.percent)),
                     Variance = pca.var.percent),
          file = "../01.1_pca/PCA_variance.csv",
          row.names = FALSE)


## -------------------------------
## 3. Set custom colors 
##    (levels must match exactly)
## -------------------------------
genotype_colors <- c("WT" = "orange", 
                     "KO" = "darkblue", 
                     "mRE" = "brown")

# Ensure factor order matches color vector
pca.dat$genotype <- factor(pca.dat$genotype, levels = names(genotype_colors))

pca_plot <- ggplot(pca.dat, aes(PC1, PC2, color = genotype)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = samples),
    size = 5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 50,
    arrow = arrow(length = unit(0.02, "npc")),
    min.segment.length = 0.1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = genotype_colors) +
  labs(
    x = paste0("PC1: ", pca.var.percent[1], "%"),
    y = paste0("PC2: ", pca.var.percent[2], "%"),
    color = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    text = element_text(family = "opensans"),
    axis.text  = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "gray95", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.text  = element_text(face = "bold")
  )




print(pca_plot)



# --- Scree Plot Data ---
cumulative_var <- round(cumsum(pca.var.percent), 2)
scree_data <- data.frame(
  PC = factor(paste0("PC", 1:length(pca.var.percent)), levels = paste0("PC", 1:length(pca.var.percent))),
  Variance = pca.var.percent,
  Cumulative = cumulative_var
)

# Scree bar plot
scree_bar_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  labs(x = NULL, y = "Variance Explained (%)") +
  theme_minimal() +
  theme(
    text = element_text(family = "opensans", size = 16),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 12.5, face = "bold"),
    axis.text.x = element_text(size = 12.5, angle = 90, hjust = 1),
    panel.background = element_rect(fill = "gray95", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )+ scale_y_continuous(
    breaks = seq(0, 25, by = 5),
    labels = function(x) paste0(x, "%"),
    limits = c(0, 25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y  = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )

# =======================================================
# OPTIONAL
# =======================================================

# Cumulative variance line plot
# scree_line_plot <- ggplot(scree_data, aes(x = PC, y = Cumulative, group = 1)) +
#   geom_line(color = "darkblue", size = 1) +
#   geom_point(color = "darkblue", size = 2) +
#   labs(x = NULL, y = "Cumulative Variance (%)") +
#   theme_minimal() +
#   theme(
#     text = element_text(family = "opensans", size = 16),
#     axis.text = element_text(size = 11, face = "bold"),
#     axis.title = element_text(size = 16, face = "bold"),
#     axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
#     panel.background = element_rect(fill = "gray95", color = NA),
#     plot.background  = element_rect(fill = "white", color = NA)
#   )

# Stack vertically
# stacked_plot <- scree_bar_plot / scree_line_plot
# print(stacked_plot)

print(scree_bar_plot)



# =======================================================
# FIGURE S1
# =======================================================


norm_exp_boxplot <- ggplot(expr_normalized_df, aes(x = Sample, y = Expression, fill = Genotype)) +
  
  # Boxplot
  geom_boxplot(alpha=0.9)  +
  
  # Same fill colors as before
  scale_fill_manual(values = c("WT" = "orange", "KO" = "darkblue", "mRE" = "brown")) +
  
  
  labs(
    title = NULL,
    x = "Samples",
    y = "Expression",
    fill = NULL
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(family = "opensans", size = 14),
    axis.text  = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )
print(norm_exp_boxplot)



# 4. Network construction by WGCNA

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices


log_mean_k <- log10(sft.data$mean.k.)

# --------------------------
# Plot 1: Scale-free topology fit
# --------------------------
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 2) +
  geom_text(nudge_y = 0.035) +
  geom_hline(yintercept = 0.8, color = 'red', linewidth = 1) +
  labs(
    x = 'Power',
    y = 'Scale-free topology model fit (signed R²)'
  ) +
  common_theme

print(a1)


# --------------------------
# Plot 2: Mean connectivity
# --------------------------
a2 <- ggplot(sft.data, aes(Power, log_mean_k, label = Power)) +
  geom_point(size = 2) +
  geom_text(nudge_y = 0.11) +
  labs(
    x = 'Power',
    y = 'Mean Connectivity (log10)'
  ) +
  common_theme

print(a2)


# Combine
combined_sft_plot <- a1 / a2
print(combined_sft_plot)



# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          networkType = "unsigned",
                          power = soft_power,
                          deepSplit = 2,
                          minModuleSize = 50,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


## ------------------------------------------------------------------
## Prepare module sizes
## ------------------------------------------------------------------
module_eigengenes <- bwnet$MEs
head(module_eigengenes)

# Number of genes per module
genes_per_module <- as.data.frame(table(bwnet$colors)) %>%
  dplyr::rename(Module = Var1, Count = Freq) %>%
  dplyr::arrange(desc(Count))

print(genes_per_module)

# Remove grey module
genes_per_module_filtered <- genes_per_module %>%
  dplyr::filter(Module != "grey")

# Colors for modules
module_colors <- setNames(
  unique(bwnet$colors[bwnet$colors != "grey"]),
  unique(bwnet$colors[bwnet$colors != "grey"])
)


# =======================================================
# FIGURE 1C
# =======================================================


## ------------------------------------------------------------------
## PLOT 1: Module size bar plot (using common_theme)
## ------------------------------------------------------------------
module_plot <- ggplot(
  genes_per_module_filtered,
  aes(x = reorder(Module, -Count), y = Count, fill = Module)
) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = module_colors) +
  scale_y_continuous(
    breaks = seq(0, max(genes_per_module_filtered$Count), by = 1000)
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 7,fontface="bold") +   # ↑ label size
  labs(
    x = "Modules",
    y = "No. of Genes"
  ) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    # GLOBAL TEXT
    text = element_text(family = "opensans", size = 22),
    
    # AXIS TITLES
    axis.title = element_text(size = 22, face = "bold"),
    
    # AXIS TICKS
    axis.text = element_text(size = 22, face = "bold"),
    
    # X AXIS LABELS
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # PANEL
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )


# Print bar plot
print(module_plot)



## ------------------------------------------------------------------
## PLOT 2: Dendrogram (theme applied where possible)
## ------------------------------------------------------------------
# Set base R plot parameters
par(
  cex.lab = 1.4,      # axis label size
  font.lab = 2,       # 2 = bold
  cex.axis = 1.2,     # axis tick size (if ticks existed)
  font.axis = 2,      # tick font bold
  cex = 1.2           # general text size inside plot
)

plotDendroAndColors(
  bwnet$dendrograms[[1]],
  bwnet$colors,
  "Module colors",
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05,
  colorHeight = 0.2,
  main = NULL
)

# Reset par to default after plot (optional)
par(mfrow = c(1,1))

# =======================================================
# FIGURE 1D
# =======================================================

# 6A. Relate modules to traits --------------------------------------------------

# binarize categorical variables

# Convert genotype to numeric indicators
colData$genotype_WT <- ifelse(colData$genotype == "WT", 1, 0)
colData$genotype_KO <- ifelse(colData$genotype == "KO", 1, 0)
colData$genotype_mRE <- ifelse(colData$genotype == "mRE", 1, 0)

colnames(colData) <- gsub("genotype_", "Genotype-", colnames(colData))


traits <- colData[, c("Genotype-WT", "Genotype-KO", "Genotype-mRE")]

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
print (heatmap.data)
head(heatmap.data)
colnames(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
heatmap.data <- heatmap.data[, colnames(heatmap.data) != "MEgrey"]

heatmap_df<- as.data.frame(heatmap.data)
print(heatmap_df)
write.csv(heatmap_df,
          file = "../08_module_trait_relationship_binarized/heatmap_data.csv",
          row.names = TRUE)

# Remove "ME" and "Genotype-" prefixes from column names
colnames(heatmap.data) <- gsub("^ME", "", colnames(heatmap.data))
colnames(heatmap.data) <- gsub("^Genotype-", "", colnames(heatmap.data))

# Check the result
print(colnames(heatmap.data))
print(MEs0$MEblack)

desired_module_order <- c(
  "black",
  "red",
  "green",
  "yellow",
  "brown",
  "blue",
  "turquoise"
)
print(colnames(heatmap.data))
trait_cols <- c("WT", "KO", "mRE")

heatmap.data <- heatmap.data[
  ,
  c(
    desired_module_order[desired_module_order %in% colnames(heatmap.data)],
    trait_cols
  ),
  drop = FALSE
]

# --- Visualization ---

# Define a single scale factor for all text elements
text_scale <- 2

p <- CorLevelPlot(
  heatmap.data,
  x = names(heatmap.data)[8:10],
  y = names(heatmap.data)[1:7],
  col = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"),
  
  # Setting all text sizes to the same scale factor
  cexTitleX = text_scale,
  cexTitleY = text_scale,
  cexLabX = text_scale,
  cexLabY = text_scale,
  cexCorval = text_scale,
  cexLabColKey = text_scale,
  
  # Setting font weights to bold (2)
  fontTitleX = 2,
  fontTitleY = 2,
  fontLabX = 2,
  fontLabY = 2,
  fontCorval = 2
)

print(p)

# =======================================================
# SUPPLEMENTARY FILE 3 - GS and MMDifferential_gene_expression.R
# =======================================================

# 6B. Intramodular analysis: Identifying driver genes ---------------

traits <- colData[, c("Genotype-WT", "Genotype-KO", "Genotype-mRE")]

# Get module names
modNames <- substring(names(module_eigengenes), 3)

for (trait in names(traits)) {
  
  # Extract the trait data
  trait_data <- as.data.frame(colData[[trait]])
  names(trait_data) <- trait
  
  # Compute module membership and significance
  geneModuleMembership <- as.data.frame(cor(norm.counts, module_eigengenes, use = 'p'))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  
  geneTraitSignificance <- as.data.frame(cor(norm.counts, trait_data, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste("GS.", trait, sep = "")
  names(GSPvalue) <- paste("p.GS.", trait, sep = "")
  
  # Combine all p-values to find the smallest non-zero value
  all_pvals <- c(as.matrix(GSPvalue), as.matrix(MMPvalue))
  min_nonzero <- min(all_pvals[all_pvals > 0])
  
  # Replace 0s with the smallest non-zero value
  GSPvalue[GSPvalue == 0] <- min_nonzero
  MMPvalue[MMPvalue == 0] <- min_nonzero
  
  
  
  # Create a dataframe with annotation and modules
  probes <- colnames(norm.counts)
  #probes2annot <- match(probes, annot$ILMN_ID)
  
  geneInfo0 <- data.frame( probes,
                           moduleColor = moduleColors,
                           geneTraitSignificance,
                           GSPvalue)
  
  # Order modules by significance for the trait
  modOrder <- order(-abs(cor(module_eigengenes, trait_data, use = 'p')))
  
  # Add module membership information
  for (mod in 1:ncol(geneModuleMembership)) {
    oldnames <- names(geneInfo0)
    geneInfo0 <- data.frame(geneInfo0,
                            geneModuleMembership[, modOrder[mod]],
                            MMPvalue[, modOrder[mod]])
    names(geneInfo0) <- c(oldnames,
                          paste("MM.", modNames[modOrder[mod]], sep = ""),
                          paste("p.MM.", modNames[modOrder[mod]], sep = ""))
  }
  
  
  
  # Rearrange data and save to CSV
  #geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0[[paste0("GS.", gsub(" ", ".", trait))]]))
  #geneInfo <- geneInfo0[geneOrder,]
  
  write.csv(geneInfo0, file = paste0("set/output/path", gsub(" ", "_", tolower(trait)), ".csv"), row.names = FALSE)
  write.table(geneInfo0, 
              file = paste0("set/output/path", 
                            gsub(" ", "_", tolower(trait)), ".tsv"), 
              row.names = FALSE, 
              sep = "\t",   # Use tab as separator
              quote = FALSE)  # Remove quotes from strings
}




selectedModule <- "turquoise"

traits <- c("Genotype-WT", "Genotype-KO", "Genotype-mRE")
traitLabels <- c("WT", "KO", "mRE")

moduleGenes <- moduleColors == selectedModule
moduleCol   <- match(selectedModule, modNames)

MM <- abs(geneModuleMembership[moduleGenes, moduleCol])

GS <- sapply(traits, function(tr)
  abs(cor(norm.counts, colData[[tr]], use = "p")[moduleGenes])
)

# ---- layout: 1 row, 3 columns ----
layout(matrix(1:3, nrow = 1))

par(
  mar = c(5, 5, 4, 2),
  oma = c(1, 1, 1, 1),
  font = 2  
)

# ---- single plotting block ----
invisible(
  mapply(
    function(gs, lab) {
      verboseScatterplot(
        x = MM,
        y = gs,
        xlab = paste("Module Membership"),
        ylab = "Gene Significance",
        main = lab,
        col = selectedModule,
        cex.main = 2.25,
        cex.lab  = 2.25,
        cex.axis = 2.25,
        xlim = c(0, 1),
        ylim = c(0, 1)
      )
    },
    as.data.frame(GS),
    traitLabels
  )
)


# 
# ---Module-wise network extraction----------
# Calculate TOM
# Step 1: Compute Topological Overlap Matrix (TOM)
print("Calculating TOM similarity matrix...")
TOM = TOMsimilarityFromExpr(norm.counts, power = 18)
print("TOM calculation completed.\n")

# Output base directory
baseDir = "set/output/path"

# Get all module colors
allModules = unique(moduleColors)
print(allModules)
print(paste("Modules detected:", paste(allModules, collapse = ", ")))

# Step 2: Calculate module eigengenes
print("Calculating module eigengenes...")
MEs = moduleEigengenes(norm.counts, moduleColors)$eigengenes
print("Module eigengenes calculated.\n")

# Step 3: Compute kME values (correlation of gene expression with module eigengenes)
print("Calculating kME (module membership) for all genes...")
kMEall = as.data.frame(cor(norm.counts, MEs, use = "p"))
print("kME calculation completed.\n")
print(moduleGenes)


# Step 4: Loop over each module
for (module in allModules) {
  
  print(paste0("Processing module: ", module, "..."))
  
  # Genes in the module
  inModule = moduleColors == module
  modProbes = colnames(norm.counts)[inModule]
  
  
  # Subset TOM
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Create module-specific output directory
  moduleDir = file.path(baseDir, module)
  if (!dir.exists(moduleDir)) {
    dir.create(moduleDir, recursive = TRUE)
    print(paste("Created directory:", moduleDir))
  }
  
  # Extract kME values for the current module
  kMEcol = paste0("ME", module)
  if (!(kMEcol %in% colnames(kMEall))) {
    warning(paste("Warning: kME column", kMEcol, "not found. Skipping module."))
    next
  }
  modkME = kMEall[modProbes, kMEcol]
  
  # Export network to Cytoscape
  print(paste("Exporting Cytoscape files for module:", module))
  exportNetworkToCytoscape(modTOM,
                           edgeFile = file.path(moduleDir, paste0("CytoscapeInput-edges-", module, ".txt")),
                           nodeFile = file.path(moduleDir, paste0("CytoscapeInput-nodes-", module, ".txt")),
                           weighted = TRUE,
                           threshold = 0.05,
                           nodeNames = modProbes,
                           altNodeNames = modProbes,
                           nodeAttr = data.frame(module = module, kME = modkME)
  )
  print(paste("Export completed for module:", module, "\n"))
}

print("Wohoooooooo! All modules processed successfully.")

## Extract the gene list per module

modules <- unique(moduleColors)

output_dir_mod_genes <- "set/output/path"

for (m in modules) {
  
  # Extract genes for the module
  genes <- colnames(norm.counts)[moduleColors == m]
  print(genes)
  
  # Create a data frame with column name
  df_genes <- data.frame(gene_ids = genes)
  
  # Save
  write.table(
    df_genes,
    file = file.path(output_dir_mod_genes, paste0("genes_", m, ".txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,   # <-- adds "gene_ids"
    sep = "\t"
  )
}
