library(dplyr)

setwd("/data/RUNX_KO_ChIPseq/")

# Different histone marks
treatment <- c("CONT", "POS")
histone_marks <- c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27ac", "H3K27me3")

for(i in 1:length(treatment)){
  for(j in 1:length(histone_marks)){
    # Import loj overlaps
    rep1_loj <- read.table(paste("loj_overlaps/", treatment[i], "_", histone_marks[j], "_rep1_loj_overlaps.bed", sep = ""), header = F, stringsAsFactors = F)
    rep2_loj <- read.table(paste("loj_overlaps/", treatment[i], "_", histone_marks[j], "_rep2_loj_overlaps.bed", sep = ""), header = F, stringsAsFactors = F)
    rep3_loj <- read.table(paste("loj_overlaps/", treatment[i], "_", histone_marks[j], "_rep3_loj_overlaps.bed", sep = ""), header = F, stringsAsFactors = F)
  
    # Import peak table
    rep1_peak_table <- read.table(paste("macs_tables/", treatment[i], "_", histone_marks[j], "_rep1_q30_peaks.xls", sep = ""), skip = 28, header = T, stringsAsFactors = F)
    rep2_peak_table <- read.table(paste("macs_tables/", treatment[i], "_", histone_marks[j], "_rep2_q30_peaks.xls", sep = ""), skip = 28, header = T, stringsAsFactors = F)
    rep3_peak_table <- read.table(paste("macs_tables/", treatment[i], "_", histone_marks[j], "_rep3_q30_peaks.xls", sep = ""), skip = 28, header = T, stringsAsFactors = F)
    
    # Generate final table 
    final_table <- matrix(ncol = 10, nrow = 1)
    colnames(final_table) <- c("chr", "start", "end", "name", "rep1_fc", "rep1_q_value", "rep2_fc", "rep2_q_value", "rep3_fc", "rep3_q_value")
  
    # Pull out unique peaks
    rep1_unique <- rep1_loj[which(rep1_loj$V9 == -1),]
    rep2_unique <- rep2_loj[which(rep2_loj$V9 == -1),]
    rep3_unique <- rep3_loj[which(rep3_loj$V9 == -1),]
    
    rep1_unique_scores <- rep1_peak_table[rep1_peak_table$name %in% rep1_unique$V4,]
    rep2_unique_scores <- rep2_peak_table[rep2_peak_table$name %in% rep2_unique$V4,]
    rep3_unique_scores <- rep3_peak_table[rep3_peak_table$name %in% rep3_unique$V4,]
    
    rep1_unique <- merge(x = rep1_unique, y = rep1_unique_scores, by.x = "V4", by.y = "name")
    rep2_unique <- merge(x = rep2_unique, y = rep2_unique_scores, by.x = "V4", by.y = "name")
    rep3_unique <- merge(x = rep3_unique, y = rep3_unique_scores, by.x = "V4", by.y = "name")
    
    if(histone_marks[j] %in% c("H3K4me2", "H3K4me3")){
      rep1_unique <- rep1_unique[,c(2,3,4,1,21,22,5,5,5,5)]
      rep2_unique <- rep2_unique[,c(2,3,4,1,5,5,21,22,5,5)]
      rep3_unique <- rep3_unique[,c(2,3,4,1,5,5,5,5,21,22)]
    } else {
      rep1_unique <- rep1_unique[,c(2,3,4,1,20,21,5,5,5,5)]
      rep2_unique <- rep2_unique[,c(2,3,4,1,5,5,20,21,5,5)]
      rep3_unique <- rep3_unique[,c(2,3,4,1,5,5,5,5,20,21)]
    }
    
    colnames(rep1_unique) <- colnames(final_table)
    colnames(rep2_unique) <- colnames(final_table)
    colnames(rep3_unique) <- colnames(final_table)
    
    final_table <- rbind(final_table, rep1_unique)
    final_table <- rbind(final_table, rep2_unique)
    final_table <- rbind(final_table, rep3_unique)
    
    final_table <- final_table[-1,]
    
    # Pull out overlapping peaks and extend coordinates 
    rep1_overlap <- rep1_loj[which(rep1_loj$V9 != -1),]
    rep2_overlap <- rep2_loj[which(rep2_loj$V9 != -1),]
    rep3_overlap <- rep3_loj[which(rep3_loj$V9 != -1),]
    
    peaks_to_add <- matrix(ncol = 10, nrow = 1)
    colnames(peaks_to_add) <- c("chr", "start", "end", "name", "rep1_fc", "rep1_q_value", "rep2_fc", "rep2_q_value", "rep3_fc", "rep3_q_value")
    
    # REPLICATE 1
    rep1_overlap_peak_names <- unique((rep1_overlap$V4))
    for(k in 1:length(rep1_overlap_peak_names)){
      message(k)
      peak_name_to_test <- rep1_overlap_peak_names[k]
      rep1_overlap_subset <- rep1_overlap[which(rep1_overlap$V4 == peak_name_to_test),]
      chromosome <- unique(rep1_overlap_subset$V1)
      start_coord <- min(c(rep1_overlap_subset$V2, rep1_overlap_subset$V9))
      end_coord <- max(c(rep1_overlap_subset$V3, rep1_overlap_subset$V10))
      peak_name_to_test_edit <- gsub(pattern = paste(treatment[i], "_", histone_marks[j], "_", sep = "") , replacement = "", x = peak_name_to_test)
      peak_overlap_names <- gsub(pattern = paste(treatment[i], "_", histone_marks[j], "_", sep = "") , replacement = "", x = rep1_overlap_subset$V11)
      peak_overlap_names <- c(peak_name_to_test_edit, peak_overlap_names)
      peak_overlap_names <- peak_overlap_names[order(peak_overlap_names)]
      peak_name <- c(paste(treatment[i], "_", histone_marks[j], "_", paste(peak_overlap_names, collapse = "_"), sep = ""))
      peak_name_split <- unlist(strsplit(peak_name, "_"))
      if(length(grep("rep1", peak_name_split)) > 0){
        rep1_idx <- grep("rep1", peak_name_split)
        rep1_score <- NULL
        rep1_q <- NULL
        for(l in 1:length(rep1_idx)){
          rep1_peak_name <- peak_name_split[rep1_idx[l]:(rep1_idx[l]+3)]
          rep1_peak_name <- paste(rep1_peak_name, collapse = "_")
          rep1_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep1_peak_name, sep = "")
          rep1_score_idx <- which(rep1_peak_table$name == rep1_peak_name)
          rep1_score <- c(rep1_score, rep1_peak_table$fold_enrichment[rep1_score_idx])
          rep1_q <- c(rep1_q, rep1_peak_table$X.log10.qvalue.[rep1_score_idx])
        }
        rep1_score <- paste(rep1_score, collapse = ";")
        rep1_q <- paste(rep1_q, collapse = ";")
      } else {
        rep1_score <- "."
        rep1_q <- "."
      }
      if(length(grep("rep2", peak_name_split)) > 0){
        rep2_idx <- grep("rep2", peak_name_split)
        rep2_score <- NULL
        rep2_q <- NULL
        for(l in 1:length(rep2_idx)){
          rep2_peak_name <- peak_name_split[rep2_idx[l]:(rep2_idx[l]+3)]
          rep2_peak_name <- paste(rep2_peak_name, collapse = "_")
          rep2_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep2_peak_name, sep = "")
          rep2_score_idx <- which(rep2_peak_table$name == rep2_peak_name)
          rep2_score <- c(rep2_score, rep2_peak_table$fold_enrichment[rep2_score_idx])
          rep2_q <- c(rep2_q, rep2_peak_table$X.log10.qvalue.[rep2_score_idx])
        }
        rep2_score <- paste(rep2_score, collapse = ";")
        rep2_q <- paste(rep2_q, collapse = ";")
      } else {
        rep2_score <- "."
        rep2_q <- "."
      }
      if(length(grep("rep3", peak_name_split)) > 0){
        rep3_idx <- grep("rep3", peak_name_split)
        rep3_score <- NULL
        rep3_q <- NULL
        for(l in 1:length(rep3_idx)){
          rep3_peak_name <- peak_name_split[rep3_idx[l]:(rep3_idx[l]+3)]
          rep3_peak_name <- paste(rep3_peak_name, collapse = "_")
          rep3_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep3_peak_name, sep = "")
          rep3_score_idx <- which(rep3_peak_table$name == rep3_peak_name)
          rep3_score <- c(rep3_score, rep3_peak_table$fold_enrichment[rep3_score_idx])
          rep3_q <- c(rep3_q, rep3_peak_table$X.log10.qvalue.[rep3_score_idx])
        }
        rep3_score <- paste(rep3_score, collapse = ";")
        rep3_q <- paste(rep3_q, collapse = ";")
      } else {
        rep3_score <- "."
        rep3_q <- "."
      }
      peaks_to_add <- rbind(peaks_to_add, c(chromosome, start_coord, end_coord, peak_name, rep1_score, rep1_q, rep2_score, rep2_q, rep3_score, rep3_q))
    }
    
    # REPLICATE 2
    rep2_overlap_peak_names <- unique((rep2_overlap$V4))
    for(k in 1:length(rep2_overlap_peak_names)){
      message(k)
      peak_name_to_test <- rep2_overlap_peak_names[k]
      rep2_overlap_subset <- rep2_overlap[which(rep2_overlap$V4 == peak_name_to_test),]
      chromosome <- unique(rep2_overlap_subset$V1)
      start_coord <- min(c(rep2_overlap_subset$V2, rep2_overlap_subset$V9))
      end_coord <- max(c(rep2_overlap_subset$V3, rep2_overlap_subset$V10))
      peak_name_to_test_edit <- gsub(pattern = paste(treatment[i], "_", histone_marks[j], "_", sep = "") , replacement = "", x = peak_name_to_test)
      peak_overlap_names <- gsub(pattern = paste(treatment[i], "_", histone_marks[j], "_", sep = "") , replacement = "", x = rep2_overlap_subset$V11)
      peak_overlap_names <- c(peak_name_to_test_edit, peak_overlap_names)
      peak_overlap_names <- peak_overlap_names[order(peak_overlap_names)]
      peak_name <- c(paste(treatment[i], "_", histone_marks[j], "_", paste(peak_overlap_names, collapse = "_"), sep = ""))
      peak_name_split <- unlist(strsplit(peak_name, "_"))
      if(length(grep("rep1", peak_name_split)) > 0){
        rep1_idx <- grep("rep1", peak_name_split)
        rep1_score <- NULL
        rep1_q <- NULL
        for(l in 1:length(rep1_idx)){
          rep1_peak_name <- peak_name_split[rep1_idx[l]:(rep1_idx[l]+3)]
          rep1_peak_name <- paste(rep1_peak_name, collapse = "_")
          rep1_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep1_peak_name, sep = "")
          rep1_score_idx <- which(rep1_peak_table$name == rep1_peak_name)
          rep1_score <- c(rep1_score, rep1_peak_table$fold_enrichment[rep1_score_idx])
          rep1_q <- c(rep1_q, rep1_peak_table$X.log10.qvalue.[rep1_score_idx])
        }
        rep1_score <- paste(rep1_score, collapse = ";")
        rep1_q <- paste(rep1_q, collapse = ";")
      } else {
        rep1_score <- "."
        rep1_q <- "."
      }
      if(length(grep("rep2", peak_name_split)) > 0){
        rep2_idx <- grep("rep2", peak_name_split)
        rep2_score <- NULL
        rep2_q <- NULL
        for(l in 1:length(rep2_idx)){
          rep2_peak_name <- peak_name_split[rep2_idx[l]:(rep2_idx[l]+3)]
          rep2_peak_name <- paste(rep2_peak_name, collapse = "_")
          rep2_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep2_peak_name, sep = "")
          rep2_score_idx <- which(rep2_peak_table$name == rep2_peak_name)
          rep2_score <- c(rep2_score, rep2_peak_table$fold_enrichment[rep2_score_idx])
          rep2_q <- c(rep2_q, rep2_peak_table$X.log10.qvalue.[rep2_score_idx])
        }
        rep2_score <- paste(rep2_score, collapse = ";")
        rep2_q <- paste(rep2_q, collapse = ";")
      } else {
        rep2_score <- "."
        rep2_q <- "."
      }
      if(length(grep("rep3", peak_name_split)) > 0){
        rep3_idx <- grep("rep3", peak_name_split)
        rep3_score <- NULL
        rep3_q <- NULL
        for(l in 1:length(rep3_idx)){
          rep3_peak_name <- peak_name_split[rep3_idx[l]:(rep3_idx[l]+3)]
          rep3_peak_name <- paste(rep3_peak_name, collapse = "_")
          rep3_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep3_peak_name, sep = "")
          rep3_score_idx <- which(rep3_peak_table$name == rep3_peak_name)
          rep3_score <- c(rep3_score, rep3_peak_table$fold_enrichment[rep3_score_idx])
          rep3_q <- c(rep3_q, rep3_peak_table$X.log10.qvalue.[rep3_score_idx])
        }
        rep3_score <- paste(rep3_score, collapse = ";")
        rep3_q <- paste(rep3_q, collapse = ";")
      } else {
        rep3_score <- "."
        rep3_q <- "."
      }
      peaks_to_add <- rbind(peaks_to_add, c(chromosome, start_coord, end_coord, peak_name, rep1_score, rep1_q, rep2_score, rep2_q, rep3_score, rep3_q))
    }
    
    # REPLICATE 3
    rep3_overlap_peak_names <- unique((rep3_overlap$V4))
    for(k in 1:length(rep3_overlap_peak_names)){
      message(k)
      peak_name_to_test <- rep3_overlap_peak_names[k]
      rep3_overlap_subset <- rep3_overlap[which(rep3_overlap$V4 == peak_name_to_test),]
      chromosome <- unique(rep3_overlap_subset$V1)
      start_coord <- min(c(rep3_overlap_subset$V2, rep3_overlap_subset$V9))
      end_coord <- max(c(rep3_overlap_subset$V3, rep3_overlap_subset$V10))
      peak_name_to_test_edit <- gsub(pattern = paste(treatment[i], "_", histone_marks[j], "_", sep = "") , replacement = "", x = peak_name_to_test)
      peak_overlap_names <- gsub(pattern = paste(treatment[i], "_", histone_marks[j], "_", sep = "") , replacement = "", x = rep3_overlap_subset$V11)
      peak_overlap_names <- c(peak_name_to_test_edit, peak_overlap_names)
      peak_overlap_names <- peak_overlap_names[order(peak_overlap_names)]
      peak_name <- c(paste(treatment[i], "_", histone_marks[j], "_", paste(peak_overlap_names, collapse = "_"), sep = ""))
      peak_name_split <- unlist(strsplit(peak_name, "_"))
      if(length(grep("rep1", peak_name_split)) > 0){
        rep1_idx <- grep("rep1", peak_name_split)
        rep1_score <- NULL
        rep1_q <- NULL
        for(l in 1:length(rep1_idx)){
          rep1_peak_name <- peak_name_split[rep1_idx[l]:(rep1_idx[l]+3)]
          rep1_peak_name <- paste(rep1_peak_name, collapse = "_")
          rep1_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep1_peak_name, sep = "")
          rep1_score_idx <- which(rep1_peak_table$name == rep1_peak_name)
          rep1_score <- c(rep1_score, rep1_peak_table$fold_enrichment[rep1_score_idx])
          rep1_q <- c(rep1_q, rep1_peak_table$X.log10.qvalue.[rep1_score_idx])
        }
        rep1_score <- paste(rep1_score, collapse = ";")
        rep1_q <- paste(rep1_q, collapse = ";")
      } else {
        rep1_score <- "."
        rep1_q <- "."
      }
      if(length(grep("rep2", peak_name_split)) > 0){
        rep2_idx <- grep("rep2", peak_name_split)
        rep2_score <- NULL
        rep2_q <- NULL
        for(l in 1:length(rep2_idx)){
          rep2_peak_name <- peak_name_split[rep2_idx[l]:(rep2_idx[l]+3)]
          rep2_peak_name <- paste(rep2_peak_name, collapse = "_")
          rep2_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep2_peak_name, sep = "")
          rep2_score_idx <- which(rep2_peak_table$name == rep2_peak_name)
          rep2_score <- c(rep2_score, rep2_peak_table$fold_enrichment[rep2_score_idx])
          rep2_q <- c(rep2_q, rep2_peak_table$X.log10.qvalue.[rep2_score_idx])
        }
        rep2_score <- paste(rep2_score, collapse = ";")
        rep2_q <- paste(rep2_q, collapse = ";")
      } else {
        rep2_score <- "."
        rep2_q <- "."
      }
      if(length(grep("rep3", peak_name_split)) > 0){
        rep3_idx <- grep("rep3", peak_name_split)
        rep3_score <- NULL
        rep3_q <- NULL
        for(l in 1:length(rep3_idx)){
          rep3_peak_name <- peak_name_split[rep3_idx[l]:(rep3_idx[l]+3)]
          rep3_peak_name <- paste(rep3_peak_name, collapse = "_")
          rep3_peak_name <- paste(treatment[i], "_", histone_marks[j], "_", rep3_peak_name, sep = "")
          rep3_score_idx <- which(rep3_peak_table$name == rep3_peak_name)
          rep3_score <- c(rep3_score, rep3_peak_table$fold_enrichment[rep3_score_idx])
          rep3_q <- c(rep3_q, rep3_peak_table$X.log10.qvalue.[rep3_score_idx])
        }
        rep3_score <- paste(rep3_score, collapse = ";")
        rep3_q <- paste(rep3_q, collapse = ";")
      } else {
        rep3_score <- "."
        rep3_q <- "."
      }
      peaks_to_add <- rbind(peaks_to_add, c(chromosome, start_coord, end_coord, peak_name, rep1_score, rep1_q, rep2_score, rep2_q, rep3_score, rep3_q))
    }
    
    peaks_to_add <- peaks_to_add[-1,]
    peaks_to_add <- as.data.frame(peaks_to_add)
    
    peaks_to_add <- distinct(peaks_to_add, chr, start, end, .keep_all = T)
    
    final_table <- rbind(final_table, peaks_to_add)
    rownames(final_table) <- NULL
    write.table(final_table, paste("replicate_overlap_results/", treatment[i], "_", histone_marks[j], "_final_overlapped_peak_table.txt", sep = ""), quote = F, sep = "\t", row.names = F)
  }
}