setwd("/RUNX1T1/RUNX_KO_ChIPseq/")

filesToEdit <- list.files(path = "consensus_overlaps/", pattern = "_overlap_table.txt", full.names = T)

histone_mark_list <- c("CONT_H3K27ac","CONT_H3K27me3","CONT_H3K4me1","CONT_H3K4me2","CONT_H3K4me3","POS_H3K27ac","POS_H3K27me3","POS_H3K4me1","POS_H3K4me2","POS_H3K4me3")

for(i in 1:length(filesToEdit)){
  message(i)
  table_import <- read.delim(file = filesToEdit[i], sep = "\t", header = T, stringsAsFactors = F)
  table_import$CONT_enhancer_status <- c(rep(NA, nrow(table_import)))
  table_import$POS_enhancer_status <- c(rep(NA, nrow(table_import)))
  for(j in 1:nrow(table_import)){
    # message(j)
    vector_of_interest <- table_import[j,11:19]
    vector_of_interest <- cbind(vector_of_interest, "1")
    colnames(vector_of_interest)[10] <- histone_mark_list[i]
    # CONT ENHANCER STATUS
    if(vector_of_interest$CONT_H3K27ac == "1" & vector_of_interest$CONT_H3K4me1 == "1" & vector_of_interest$CONT_H3K27me3 == "0" & vector_of_interest$CONT_H3K4me3 == "0"){
      table_import$CONT_enhancer_status[j] <- "Active"
    } else if(vector_of_interest$CONT_H3K27me3 == "1" & vector_of_interest$CONT_H3K4me1 == "1" & vector_of_interest$CONT_H3K27ac == "0" & vector_of_interest$CONT_H3K4me3 == "0"){
      table_import$CONT_enhancer_status[j] <- "Poised"
    } else if(vector_of_interest$CONT_H3K4me1 == "1" & vector_of_interest$CONT_H3K27ac == "0" & vector_of_interest$CONT_H3K27me3 == "0" & vector_of_interest$CONT_H3K4me3 == "0"){
      table_import$CONT_enhancer_status[j] <- "Primed"
    } else {
      table_import$CONT_enhancer_status[j] <- "None/Insufficient information"
    }
    
    # POS ENHANCER STATUS
    if(vector_of_interest$POS_H3K27ac == "1" & vector_of_interest$POS_H3K4me1 == "1" & vector_of_interest$POS_H3K27me3 == "0" & vector_of_interest$POS_H3K4me3 == "0"){
      table_import$POS_enhancer_status[j] <- "Active"
    } else if(vector_of_interest$POS_H3K27me3 == "1" & vector_of_interest$POS_H3K4me1 == "1" & vector_of_interest$POS_H3K27ac == "0" & vector_of_interest$POS_H3K4me3 == "0"){
      table_import$POS_enhancer_status[j] <- "Poised"
    } else if(vector_of_interest$POS_H3K4me1 == "1" & vector_of_interest$POS_H3K27ac == "0" & vector_of_interest$POS_H3K27me3 == "0" & vector_of_interest$POS_H3K4me3 == "0"){
      table_import$POS_enhancer_status[j] <- "Primed"
    } else {
      table_import$POS_enhancer_status[j] <- "None/Insufficient information"
    }
    
  }
  write.table(table_import, paste("consensus_overlaps/", histone_mark_list[i],"_annotated_enhancer_status_overlap_table.txt", sep = ""), sep = "\t", quote = F, row.names = F)
}
