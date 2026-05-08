############################################################
# IGF2BP2 PROJECT – WGCNA
############################################################

# Author: Vishnu A M
# Lab: CLOCK Lab
# Institute: Indian Institute of Technology, Gandhinagar, India
#
# This script runs WGCNA on the Tximport generated raw gene expression matrix

# Input:
# Processed BigWig files
#
# ### Output Figures #
# ChIP-seq plot for the region of interest

############################################################
# 1️⃣ Working Directory
############################################################
setwd("path/to/working/directory")

############################################################
# 2️⃣ Libraries
############################################################
library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)

############################################################
# 3️⃣ Genomic Region
############################################################
chr  <- "chr16"
from <- 22160000
to   <- 22165000

############################################################
# 4️⃣ Your Four Colors
############################################################
cols <- list(
  "ChIP ZT10"  = "#185b8c",
  "ChIP ZT22"  = "#73b7e7",
  "Input ZT10" = "#d16200",
  "Input ZT22" = "#ffa85c"
)

y_limits <- c(0, 60)

############################################################
# 5️⃣ Genome Axis (Converted to Kbs)
############################################################
gtrack <- GenomeAxisTrack(
  col = "black",
  fontcolor = "black",
  cex = 0.75,                  # Increased from 0.35
  labelPos = "alternating",   # Prevents overlap
  exponent = 3,               # Divides ticks by 1,000 to show as Kb
  name = "Location (Kb)",     # Labels the track
  showTitle = TRUE
)

############################################################
# 6️⃣ BigWig Tracks
############################################################
tr_zt10 <- DataTrack(range="TreatmentZT10.bw", type="histogram", genome="mm10", chromosome=chr, name="ChIP ZT10", fill.histogram=cols[["ChIP ZT10"]], col.histogram=cols[["ChIP ZT10"]], ylim=y_limits)
tr_zt22 <- DataTrack(range="TreatmentZT22.bw", type="histogram", genome="mm10", chromosome=chr, name="ChIP ZT22", fill.histogram=cols[["ChIP ZT22"]], col.histogram=cols[["ChIP ZT22"]], ylim=y_limits)
in_zt10 <- DataTrack(range="InputZT10.bw", type="histogram", genome="mm10", chromosome=chr, name="Input ZT10", fill.histogram=cols[["Input ZT10"]], col.histogram=cols[["Input ZT10"]], ylim=y_limits)
in_zt22 <- DataTrack(range="InputZT22.bw", type="histogram", genome="mm10", chromosome=chr, name="Input ZT22", fill.histogram=cols[["Input ZT22"]], col.histogram=cols[["Input ZT22"]], ylim=y_limits)

# Limit Y-ticks to only 0 and 60 so they take up less space
displayPars(tr_zt10) <- list(showAxis = TRUE, yTicksAt = c(0, 60))
displayPars(tr_zt22) <- list(showAxis = TRUE, yTicksAt = c(0, 60))
displayPars(in_zt10) <- list(showAxis = TRUE, yTicksAt = c(0, 60))
displayPars(in_zt22) <- list(showAxis = TRUE, yTicksAt = c(0, 60))

# ############################################################
# # 7️⃣ Minimal ENCODE Block
# ############################################################
# encode_gr <- GRanges(seqnames=chr, ranges=IRanges(start=22163339, end=22163539))
# 
# encode_track <- AnnotationTrack(
#   encode_gr, genome="mm10", chromosome=chr,
#   name="ENCODE cCREs", fill="red", col="red", shape="box"
# )
# 
# displayPars(encode_track) <- list(
#   cex.title = 0.5, # Increased from 0.25
#   cex = 0.5        # Increased from 0.25
# )

############################################################
# 8️⃣ Combine Tracks (RefSeq track removed!)
############################################################
track_list <- list(gtrack, tr_zt10, tr_zt22, in_zt10, in_zt22)

# Updated sizes vector to match 6 tracks instead of 7
sizes_clean <- c(0.6, 1.4, 1.4, 1.4, 1.4)

############################################################
# 9️⃣ SAFE PREVIEW (PDF)
############################################################
pdf("preview.pdf", width=5, height=4)
plotTracks(
  track_list, from=from, to=to, sizes=sizes_clean,
  background.panel="white", background.title="white",
  col.title="black", col.axis="black", col.border.title=NA,
  rotation.title=90,
  title.width=0.9,    # Slightly widened to accommodate larger font
  innerMargin=2,      # SQUEEZES the tracks closer together
  cex.title=0.6,      # Main BigWig track titles increased from 0.4
  cex.axis=0.45       # Y-axis numbers increased from 0.3
)
dev.off()

############################################################
# 🔟 Final 300 DPI TIFF (Manuscript)
############################################################
tiff("nr1d1_igf2bp2_chip.tiff",
     width=5, height=5, units="in", res=300, compression="lzw", bg="transparent")

png("nr1d1_igf2bp2_chip.png",
     width=5, height=5, units="in", res=300, bg="transparent")

plotTracks(
  track_list, from=from, to=to, sizes=sizes_clean,
  background.panel="transparent", background.title="transparent",
  col.title="black", col.axis="black", col.border.title=NA,
  rotation.title=90,
  title.width=0.6,    
  innerMargin=1,      
  cex.title=0.85,      
  cex.axis=0.5        
)
dev.off()
