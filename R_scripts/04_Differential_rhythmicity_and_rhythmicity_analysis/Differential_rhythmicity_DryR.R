# ============================================================
# IGF2BP2 PROJECT – DryR Circadian Analysis
# ============================================================

# Author: Vishnu A M
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script performs Differential rhythmicity analysis using DryR

# Input:
# DDS75 filtered raw gene expression matrix
#
# Output:
# DryR models with the differential rhythimcity parameters for the all the input genes.
# Figure 2B was plotted using the output of this script in python.
#



# ============================================================
# IGF2BP2 PROJECT – DryR Circadian Analysis
# Clean & Reproducible Version
# ============================================================

#save the workspace image
save.image("path/to//.RData")


# 1. Load the data
#load the workspace image
load("/path/to/.RData")



# -------------------------------
# 1. Load libraries
# -------------------------------
library(tidyverse)
library(dryR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpattern)

# -------------------------------
# 2. Load raw count data
# -------------------------------
countData <- read.delim(
  "path/to//dds75_filtered_raw_counts.tsv",
  stringsAsFactors = FALSE
)

# Set Ensembl IDs as rownames
rownames(countData) <- countData$X
countData$X <- NULL

# Clean column names
colnames(countData) <- gsub("\\.", "_", colnames(countData))

# -------------------------------
# 3. Map Ensembl → Gene Symbol
# -------------------------------
gene_map <- bitr(
  rownames(countData),
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Mm.eg.db
)

# Remove duplicate ENSEMBL mappings
gene_map_unique <- gene_map %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# -------------------------------
# 4. Merge counts + gene symbols
# -------------------------------
countData_merged <- countData %>%
  mutate(ENSEMBL = rownames(countData)) %>%
  left_join(gene_map_unique, by = "ENSEMBL") %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", ENSEMBL, SYMBOL))

# Define expression columns explicitly
expr_cols <- colnames(countData)

# Keep highest expressed gene if duplicate SYMBOL exists
countData_final <- countData_merged %>%
  rowwise() %>%
  mutate(total_expr = sum(c_across(all_of(expr_cols)))) %>%
  ungroup() %>%
  arrange(desc(total_expr)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(-total_expr)

# -------------------------------
# 5. Prepare matrix
# -------------------------------
countData_final <- as.data.frame(countData_final)

rownames(countData_final) <- countData_final$SYMBOL
countData_final$SYMBOL  <- NULL
countData_final$ENSEMBL <- NULL

# Round counts
countData_final <- round(countData_final)

# Dynamic filtering (≥15 counts in ≥75% samples)
keep <- rowSums(countData_final >= 15) >= 0.75 * ncol(countData_final)
countData_dryr <- countData_final[keep, ]

# Manual fix if needed
rownames(countData_dryr)[rownames(countData_dryr) == "ENSMUSG00000113216"] <- "Gm40841"

# -------------------------------
# 6. Define experimental design
# -------------------------------
time_points <- c(0, 4, 8, 12, 16, 20)
time_vector <- rep(time_points, each = 3)
time <- rep(time_vector, times = 3)

group <- sub("_ZT.*", "", colnames(countData_dryr))
group <- factor(group)

# -------------------------------
# 7. Run DryR
# -------------------------------
dryList <- dryseq(countData_dryr, group, time)

dry_res <- dryList$results %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene")

# Save results
write.table(
  dry_res,
  file = "path/to/dryr_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 8. Rhythmicity Classification
# ============================================================

classify_condition <- function(mean_ref, amp_ref, phase_ref,
                               mean_x,   amp_x,   phase_x) {
  
  if (any(is.na(c(mean_ref, amp_ref, phase_ref,
                  mean_x,   amp_x,   phase_x)))) {
    return("arrhythmic")
  }
  
  amp_same   <- abs(amp_x - amp_ref) < 1e-6
  phase_same <- abs(phase_x - phase_ref) < 1e-6
  
  if (amp_same & phase_same) {
    "normal"
  } else if (!amp_same & phase_same) {
    "alt_amp"
  } else if (amp_same & !phase_same) {
    "alt_phase"
  } else {
    "alt_phase_amp"
  }
}

status_df <- dry_res %>%
  rowwise() %>%
  mutate(
    WT = ifelse(
      any(is.na(c(mean_WT, amp_WT, phase_WT))),
      "arrhythmic",
      "normal"
    ),
    KO = classify_condition(
      mean_WT, amp_WT, phase_WT,
      mean_KO, amp_KO, phase_KO
    ),
    mRE = classify_condition(
      mean_WT, amp_WT, phase_WT,
      mean_mRE, amp_mRE, phase_mRE
    )
  ) %>%
  ungroup() %>%
  dplyr::select(Gene, WT, KO, mRE)

# ============================================================
# 9. Hub Gene Analysis
# ============================================================

# --- Consistent Font Size ---
font_size_global <- 20 

# 1. Clean data and define explicit factor levels
plot_df <- status_hubs_wt %>%
  pivot_longer(
    cols = c(WT, KO, mRE),
    names_to = "Condition",
    values_to = "Rhythmicity"
  ) %>%
  mutate(
    Rhythmicity = case_when(
      Rhythmicity %in% c("alt_phase", "alt_amp") ~ "alt_phase_amp",
      TRUE ~ Rhythmicity
    ),
    Rhythmicity = factor(Rhythmicity, levels = c("normal", "alt_phase_amp", "arrhythmic")),
    Condition = factor(Condition, levels = c("WT", "KO", "mRE"))
  )

# 2. Define the Mappings (Solid colors only)
fill_colors <- c(
  "normal" = "orange",
  "alt_phase_amp" = "blue",
  "arrhythmic" = "white"
)

############################################################
# This output file was then used for plotting FIGURE 2B in python
############################################################


# save the plotting file

status_hubs_wt_df<- as.data.frame(status_hubs_wt)

write.table(
  status_hubs_wt_df,
  file = "path/to//dryr_hubs_rhytmicity.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 3. Plot
p <- ggplot(
  plot_df,
  aes(
    x = Condition,
    y = Gene,
    fill = Rhythmicity
  )
) +
  # Using standard geom_tile for a clean, solid-cell look
  geom_tile(
    color = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = fill_colors,
    breaks = names(fill_colors),
    labels = c("normal" = "Normal", "alt_phase_amp" = "Altered", "arrhythmic" = "Arrhythmic")
  ) +
  coord_fixed() +
  theme_minimal(base_size = font_size_global) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      angle = 45,
      hjust = 1,
      colour = "black",
      size = font_size_global
    ),
    axis.text.y = element_text(
      face = "bold",
      colour = "black",
      size = font_size_global
    ),
    legend.text = element_text(face = "bold", size = font_size_global),
    legend.title = element_text(face = "bold", size = font_size_global),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(fill = "Rhythmicity")

print(p)
