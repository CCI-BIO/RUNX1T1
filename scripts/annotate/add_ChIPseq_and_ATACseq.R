library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

# Set working directory 
setwd(dir = "/RUNX1T1/RUNX_KO_ChIPseq/")

# List of files
fileList <- list.files(path = "consensus_overlaps/", pattern = "annotated_enhancer_status_overlap_table.txt", full.names = T)
fileNames <- list.files(path = "consensus_overlaps/", pattern = "annotated_enhancer_status_overlap_table.txt", full.names = F)
fileNames <- gsub(pattern = ".txt", replacement = "", x = fileNames)

repFiles <- list.files(path = "Giorgio_ChIPseq_plus_ATACseq_bed/", pattern = ".bed", full.names = T)

# Determine overlaps using Genomic Ranges 
for(i in 1:length(fileList)){
  previous_table <- read.delim(fileList[i], header = T, stringsAsFactors = F)
  previous_coords <- previous_table[,c(1,2,3)]
  previous_coords_granges <- makeGRangesFromDataFrame(df = previous_coords, seqnames.field = "chrom")
  previous_table$LSD1SH1_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$LSD1SH2_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$LSD1SH3_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$LSD1SH4_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$RCOR3_2_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$RCOR3_3_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$RCOR3_A_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$RUNX1T1_rep1_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$RUNX1T1_rep2_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$RUNX1T1_rep3_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$shAMB_1_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$shAMB_2_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$shAMB_Doxy_1_peaks <- c(rep(0, nrow(previous_table)))
  previous_table$shAMB_Doxy_2_peaks <- c(rep(0, nrow(previous_table)))
  for(j in 1:length(repFiles)){
    message(j)
    binding_table <- read.delim(repFiles[j], header = F, stringsAsFactors = F)
    binding_coords <- binding_table[,c(1,2,3)]
    colnames(binding_coords) <- c("chrom", "start", "end")
    binding_coords_granges <- makeGRangesFromDataFrame(df = binding_coords, seqnames.field = "chrom")
    overlap_regions <- findOverlaps(query = binding_coords_granges, subject = previous_coords_granges)
    previous_table[subjectHits(overlap_regions),(21+j)] <- 1
  }
  write.table(x = previous_table, file = paste("consensus_overlaps/", fileNames[i], "_added_ChIPseq_and_ATACseq.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}
