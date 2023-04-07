library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

setwd("/RUNX1T1/RUNX_KO_ChIPseq/")

symbol <- keys(org.Hs.eg.db, "SYMBOL")
g <- select(org.Hs.eg.db, columns = c("SYMBOL", "ENTREZID"), keys = symbol, keytype = "SYMBOL")
chrloc <- select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys = g$ENTREZID, keytype = "GENEID", columns = c("TXCHROM", "TXSTART", "TXEND"))
for(i in 1:nrow(chrloc)){
  chrloc$GENENAME[i] <- g$SYMBOL[which(g$ENTREZID == chrloc$GENEID[i])]
}

# Import all peaks
fileList <- list.files("consensus_overlaps/", pattern = "overlaps.bed", full.names = T)

summitFileList <- list.files("replicate_overlap_results/", pattern = ".bed", full.names = T)

all_treatment_labels <- c("CONT_H3K27ac","CONT_H3K27me3","CONT_H3K4me1","CONT_H3K4me2","CONT_H3K4me3","POS_H3K27ac","POS_H3K27me3","POS_H3K4me1","POS_H3K4me2","POS_H3K4me3")

for(i in 1:length(fileList)){
  message(i)
  peak_table <- read.table(file = summitFileList[i], header = F, sep = "\t", stringsAsFactors = F)
  table_import <- read.table(file = fileList[i], header = F, sep = "\t", stringsAsFactors = F)
  final_table <- matrix(nrow = nrow(peak_table), ncol = 19)
  colnames(final_table) <- c("chrom","start","end","id","score","gene_5k","gene_50k","gene_100k","gene_200k","gene_1mb",all_treatment_labels[-i])
  final_table <- as.data.frame(final_table)
  final_table[,c(1:5)] <- peak_table[,c(1,2,3,4,5)]
  for(j in 1:nrow(final_table)){
    idToTest <- final_table$id[j]
    table_rows <- table_import[which(table_import[,4] == idToTest),]
    if(nrow(table_rows) > 0){
      vec <- c(rep(0, 9))
      idx <- which(all_treatment_labels[-i] %in% table_rows[,7])
      vec[idx] <- 1
    } else {
      vec <- c(rep(0, 9))
    }
    final_table[j,11:ncol(final_table)] <- vec
    chromToTest <- final_table$chrom[j]
    # 5k
    start_5k <- final_table$start[j] - 5000
    end_5k <- final_table$end[j] + 5000
    subset_5k <- chrloc[which(chrloc$TXCHROM == chromToTest & start_5k <= chrloc$TXSTART & chrloc$TXSTART <= end_5k),]
    if(nrow(subset_5k) > 0){
      final_table$gene_5k[j] <- paste(unique(subset_5k$GENENAME), collapse = " | ")
    } else if(nrow(subset_5k) == 0){
      final_table$gene_5k[j] <- NA
    }
    # 50k
    start_50k <- final_table$start[j] - 50000
    end_50k <- final_table$end[j] + 50000
    subset_50k <- chrloc[which(chrloc$TXCHROM == chromToTest & start_50k <= chrloc$TXSTART & chrloc$TXSTART <= end_50k),]
    if(nrow(subset_50k) > 0){
      final_table$gene_50k[j] <- paste(unique(subset_50k$GENENAME), collapse = " | ")
    } else if(nrow(subset_50k) == 0){
      final_table$gene_50k[j] <- NA
    }
    # 100k
    start_100k <- final_table$start[j] - 100000
    end_100k <- final_table$end[j] + 100000
    subset_100k <- chrloc[which(chrloc$TXCHROM == chromToTest & start_100k <= chrloc$TXSTART & chrloc$TXSTART <= end_100k),]
    if(nrow(subset_100k) > 0){
      final_table$gene_100k[j] <- paste(unique(subset_100k$GENENAME), collapse = " | ")
    } else if(nrow(subset_100k) == 0){
      final_table$gene_100k[j] <- NA
    }
    # 200k
    start_200k <- final_table$start[j] - 200000
    end_200k <- final_table$end[j] + 200000
    subset_200k <- chrloc[which(chrloc$TXCHROM == chromToTest & start_200k <= chrloc$TXSTART & chrloc$TXSTART <= end_200k),]
    if(nrow(subset_200k) > 0){
      final_table$gene_200k[j] <- paste(unique(subset_200k$GENENAME), collapse = " | ")
    } else if(nrow(subset_200k) == 0){
      final_table$gene_200k[j] <- NA
    }
    # 1mb
    start_1mb <- final_table$start[j] - 1000000
    end_1mb <- final_table$end[j] + 1000000
    subset_1mb <- chrloc[which(chrloc$TXCHROM == chromToTest & start_1mb <= chrloc$TXSTART & chrloc$TXSTART <= end_1mb),]
    if(nrow(subset_1mb) > 0){
      final_table$gene_1mb[j] <- paste(unique(subset_1mb$GENENAME), collapse = " | ")
    } else if(nrow(subset_1mb) == 0){
      final_table$gene_1mb[j] <- NA
    }
  }
  write.table(x = final_table, file = paste("consensus_overlaps/", all_treatment_labels[i], "_overlap_table.txt", sep = ""), sep = "\t", row.names = F, quote = F)
}
