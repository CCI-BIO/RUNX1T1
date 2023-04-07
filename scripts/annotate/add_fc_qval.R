# Add fc and qval to tables

setwd("/RUNX1T1/RUNX_KO_ChIPseq/")

finalTableList <- list.files(path = "consensus_overlaps/", pattern = "added_LSD1_RCOR_rep.txt", full.names = T)
finalTableListNames <- list.files(path = "consensus_overlaps/", pattern = "added_LSD1_RCOR_rep.txt", full.names = F)
finalTableListNames <- gsub(pattern = ".txt", replacement = "", x = finalTableListNames)
scoreBedFiles <- list.files(path = "replicate_overlap_results/", pattern = ".txt", full.names = T)

for(i in 1:length(finalTableList)){
  final_table <- read.delim(finalTableList[i], header = T, stringsAsFactors = F)
  bed_file <- read.delim(scoreBedFiles[i], header = T, stringsAsFactors = F)
  final_table <- merge(x = final_table, y = bed_file, by.x = c("chrom", "start", "end", "id"), by.y = c("chr", "start", "end", "name"), all.x = T)
  write.table(final_table, paste("consensus_overlaps/", finalTableListNames[i], "_added_scores.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}