# ============================================================
# IGF2BP2 PROJECT – Gene ontology
# ============================================================

# Author: Shriyansh Srivastava
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script performs Gene ontology analysis using the Clusterprofiler package

# Input:
# Output from the RAIN and Metacycle
#
# Output:
# Heatmap showing the rhythmicity of Igf2bp2 in other circadian datasets using RAIN and Metacycle
#






setwd('C:/Users/l/Downloads')
getwd()


library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# File path
file_path <- "submission_table_for_rhythmicity-details.xlsx"

# Get all sheet names
sheet_names <- excel_sheets(file_path)

# Read each sheet into a list of data frames
df_list <- lapply(sheet_names, function(x) read_excel(file_path, sheet = x,skip=2))

# Assign sheet names to the list
names(df_list) <- sheet_names

# Optional: create separate data frames in the environment
list2env(df_list, envir = .GlobalEnv)


metacycle <- na.omit(metacycle)
rain <- na.omit(rain)

metacycle[,4:ncol(metacycle)] <- round(metacycle[,4:ncol(metacycle)], 3)
rain[,4:ncol(rain)] <- round(rain[,4:ncol(rain)], 3)


metacycle <- metacycle %>% filter(GeneSymbol == "Igf2bp2")
rain <- rain %>% filter(GeneSymbol == "Igf2bp2")

metacycle <- metacycle %>%
  filter(GeneSymbol == "Igf2bp2") %>%
  mutate(
    Rhythmic = ifelse(meta2d_pvalue < 0.05 & meta2d_BH.Q < 0.1, #change meta2d_BH.Q to 0.05 if needed
                      "Rhythmic", "Arrhythmic"),
    Species = ifelse(grepl("GSE182117|GSE109825|GSE108539",
                           `GEO ID and Experimental group`),
                     "Human", "Mouse")
  )

rain <- rain %>%
  filter(GeneSymbol == "Igf2bp2") %>%
  mutate(
    Rhythmic = ifelse(pvalue < 0.05 & FDR < 0.1, #change FDR to 0.05 fi needed
                      "Rhythmic", "Arrhythmic"),
    Species = ifelse(grepl("GSE182117|GSE109825|GSE108539",
                           `GEO ID and Experimental group`),
                     "Human", "Mouse")
  )




# prepare metacycle
meta_plot <- metacycle %>%
  select(Dataset = `GEO ID and Experimental group`,
         Species,
         p = meta2d_pvalue,
         q = meta2d_BH.Q,
         Rhythmic) %>%
  mutate(
    Method = "MetaCycle",
    label = paste0("bold('p=", round(p,3), ", q=", round(q,3), "')")
  )

# prepare rain
rain_plot <- rain %>%
  select(Dataset = `GEO ID and Experimental group`,
         Species,
         p = pvalue,
         q = FDR,
         Rhythmic) %>%
  mutate(
    Method = "RAIN",
    label = paste0("bold('p=", round(p,3), ", q=", round(q,3), "')")
  )

# combine
df_plot <- bind_rows(meta_plot, rain_plot)


dataset_order <- unique(df_plot$Dataset)

# remove underscores
df_plot$Dataset <- gsub("_"," ", df_plot$Dataset)

# restore order and reverse (so first dataset appears at top)
df_plot$Dataset <- factor(
  df_plot$Dataset,
  levels = rev(gsub("_"," ", dataset_order))
)

# method order
df_plot$Method <- factor(df_plot$Method, levels=c("MetaCycle","RAIN"))

# species order
df_plot$Species <- factor(df_plot$Species, levels=c("Mouse","Human"))


# plot
ggplot(df_plot, aes(Method, Dataset, fill = Rhythmic)) +
  
  geom_tile(color = "grey90", linewidth = 0.6) +
  
  geom_text(aes(label = label), size = 3,parse=TRUE) +
  
  scale_fill_manual(values = c(
    "Arrhythmic" = "#E6E6E6",
    "Rhythmic" = "#1F77B4"
  )) +
  
  facet_grid(Species ~ ., scales = "free_y", space = "free_y") +
  
  labs(
    title = "Igf2bp2 Rhythmicity Across Datasets",
    x = "",
    y = "",
    fill = ""
  ) +
  
  guides(fill = guide_legend(nrow = 1)) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text.y = element_text(angle = 0, face = "bold"),
    plot.title = element_text(hjust = 1),
    legend.position = "bottom"
  ) #550 x750

# > save.image("~/files_shared_by_vishnu/igf2bp2_paper/rhythmicity_analyses/heatmapcode_from_file.RData")
