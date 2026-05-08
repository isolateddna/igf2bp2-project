
# ============================================================
# IGF2BP2 PROJECT – Gene ontology
# ============================================================

# Author: Shriyansh Srivastava
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script performs Gene ontology analysis using the Clusterprofiler package

# Input:
# Common genes between Nr1d1 and Igf2bp2 were parsed into the enrichGO fucntion
#
# Output:
# Biological pathways enriched with the enriched genes
#



# ============================================================
# FIGURE S3
# ============================================================


setwd("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/network_analysis/turquoise/Nr1d2_Igf2bp2_common_neighbors/")
getwd()
# ==============================
# GO ENRICHMENT PIPELINE
# clusterProfiler workflow
# ==============================

# 1 Install packages (run once)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("clusterProfiler","org.Mm.eg.db","enrichplot","ggplot2"))

# 2 Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

# ==============================
# 3 Gene list (your genes)
# ==============================

genes <- c(
  "Gm40841","Tcap","Per3","Clock","5930430L01Rik","Tef","Ptpn4","Fam110b",
  "Zdbf2","Coq10b","Dbp","Slc46a3","Ttc4","Mid1ip1","Wsb1","Sesn3","Gpam",
  "Wls","Slc25a33","Usp2","Nlrc3","Dusp14","4930578G10Rik","Chmp1b2","Mboat1",
  "Stk35","Slc7a2","Parp8","Nfil3","Pnpla3","Gm21955","Ypel2","Tmc7","Hcfc1r1",
  "A530013C23Rik","Cdc14b","Zfp612","Tet2","Slc25a36","Mkrn1","Bcar3","Leo1",
  "Dyrk2","Gm64932","Cry1","Amer1","Lrrc30","Asph","Ppp1r27","Miga2","Per2",
  "Sdc4","Elf2","Marchf4","Serpine1"#,"Igf2bp2","Nr1d2"
)

# ==============================
# 4 Convert SYMBOL → ENTREZ
# ==============================

gene_df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

gene_ids <- gene_df$ENTREZID

# ==============================
# 5 OPTIONAL: define background
# ==============================

# Example: if you have all expressed genes
background_genes <- read.table("../../../../output_text_files/03_module_info/turquoise_genes.txt")

background_genes <- background_genes$V1
# Convert background genes to Entrez
bg_df <- bitr(background_genes,
               fromType="ENSEMBL",
               toType="ENTREZID",
               OrgDb=org.Mm.eg.db)

bg_ids <- bg_df$ENTREZID

# If you DON'T have background, comment universe line below

# ==============================
# 6 Run GO enrichment
# ==============================

ego_bp <- enrichGO(
  gene          = gene_ids,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE,
  # universe     = bg_ids
)

# ==============================
# 7 Remove redundant GO terms
# ==============================

ego_bp_simple <- simplify(
  ego_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# ==============================
# 8 Plot results
# ==============================

dotplot(ego_bp_simple, showCategory = 30) #610 x 650

# Barplot
barplot(
  ego_bp_simple,
  showCategory = 30,
  x = "Count"
)

# GO network
cnetplot(ego_bp_simple, showCategory = 30) #1100 x 750

cnetplot(
  ego_bp_simple,
  showCategory = 10,
  node_label = "all",
  layout = "fr"
) #570 x 490


save.image(file = "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/rstudio_codes_final/go/common_neighbors_and_circadian_community_GO.RData")

# GO similarity network
# ego_bp_simple <- pairwise_termsim(ego_bp_simple)
# emapplot(ego_bp_simple)






# ==============================
# 9 Save results
# ==============================
# 
# write.csv(as.data.frame(ego_bp_simple),
#           "GO_BP_simplified_results_for_commonn_neighbors_in_circadian_community.csv")
# 
# gobp <- read.csv("GO_BP_simplified_results_for_commonn_neighbors_in_circadian_community.csv")
# # Save figure
# pdf("GO_dotplot.pdf", width=8, height=6)
# dotplot(ego_bp_simple, showCategory = 15)
# dev.off()