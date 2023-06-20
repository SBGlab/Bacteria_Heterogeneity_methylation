setwd("E:/Universidad/MBB/TFM/methylation_motifs")
library(Biostrings)
library(tidyverse)


motifs=function(motif_22_filtered,motif_22_full,methylation_full_fasta){
    #FULL DATA
  motif_full=read.csv(file=motif_22_full,header = FALSE, sep = "\t",comment=";")
  methylation_number_full=length(table(motif_full$V4))
  motif_full=motif_full[!(motif_full$V5==-22 & motif_full$V6==-1),]
  #allow a slip of 1 bp call methylation
  motif_full_filtered <- motif_full[motif_full$V1 == "GATC" & motif_full$V5 == -13 & motif_full$V6 == -10 | motif_full$V1 == "GATC" & motif_full$V5 == -14 & motif_full$V6 == -11 | motif_full$V1 != "GATC", ]
  motif_full_filtered <- motif_full_filtered[motif_full_filtered$V1 == "CCWGG" & motif_full_filtered$V5 == -13 & motif_full_filtered$V6 == -9 | motif_full_filtered$V1 == "CCWGG" & motif_full_filtered$V5 == -15 & motif_full_filtered$V6 == -11 | motif_full_filtered$V1 != "CCWGG", ]
  motif_full_filtered <- motif_full_filtered[motif_full_filtered$V1 == "AACNNNNNNGTGC" & motif_full_filtered$V5 == -13 & motif_full_filtered$V6 == -1 | motif_full_filtered$V1 == "AACNNNNNNGTGC" & motif_full_filtered$V5 == -21 & motif_full_filtered$V6 == -9 | motif_full_filtered$V1 != "AACNNNNNNGTGC", ]
  motif_full_filtered <- motif_full_filtered[motif_full_filtered$V1 == "GCACNNNNNNGTT" & motif_full_filtered$V5 == -14 & motif_full_filtered$V6 == -2 | motif_full_filtered$V1 == "GCACNNNNNNGTT" & motif_full_filtered$V5 == -23 & motif_full_filtered$V6 == -11 | motif_full_filtered$V1 != "GCACNNNNNNGTT", ]
  motif_full_filtered <- motif_full_filtered[motif_full_filtered$V1 == "ATGCAT" & motif_full_filtered$V5 == -16 & motif_full_filtered$V6 == -11 | motif_full_filtered$V1 == "ATGCAT" & motif_full_filtered$V5 == -13 & motif_full_filtered$V6 == -8 | motif_full_filtered$V1 != "ATGCAT", ]
  #list of pairs of indices that are duplicated
  duplicated_indices=which(duplicated(motif_full_filtered$V4))
  result_list <- list()
  for (num in duplicated_indices) {
    result_vector <- c(num, num - 1)
    result_list <- append(result_list, list(result_vector))
  }
  removal_indices <- c()
  # Iterate through the paired indices and compare scores in V8
  for (pair in result_list) {
    index1 <- pair[1]
    index2 <- pair[2]
    score1 <- as.numeric(motif_full_filtered[index1, "V8"])
    score2 <- as.numeric(motif_full_filtered[index2, "V8"])
    # Identify the index of the row with the highest score
    if (score1 > score2) {
      removal_indices <- c(removal_indices, index2)
    } else {
      removal_indices <- c(removal_indices, index1)
    }
  }
  # Remove the rows with the indices from the removal list
  data_full <- motif_full_filtered[-removal_indices, ]
  motif_number_full=dim(data_full)[1]
  percentage_motifs_full=motif_number_full/methylation_number_full*100
  percentage_motifs_GATC_full=sum(data_full$V1=="GATC")/motif_number_full*100
  percentage_motifs_AACNNNNNNGTGC_full=(sum(data_full$V1=="AACNNNNNNGTGC")+sum(data_full$V1=="GCACNNNNNNGTT"))/motif_number_full*100
  percentage_motifs_CCWGG_full=sum(data_full$V1=="CCWGG")/motif_number_full*100
  percentage_motifs_ATGCAT_full=sum(data_full$V1=="ATGCAT")/motif_number_full*100
  score_GATC_full=mean(unlist(filter(data_full, V1 == "GATC")%>%select(V8)))
  score_AACNNNNNNGTGC_full=mean(unlist(filter(data_full, V1 == "AACNNNNNNGTGC" | V1 == "GCACNNNNNNGTT")%>%select(V8)))
  score_CCWGG_full=mean(unlist(filter(data_full, V1 == "CCWGG")%>%select(V8)))
  score_ATGCAT_full=mean(unlist(filter(data_full, V1 == "ATGCAT")%>%select(V8)))
  
  
  #NOW FOR THE FILTERED DATA
  motif_filtered=read.csv(file=motif_22_filtered,header = FALSE, sep = "\t",comment=";")
  methylation_number_filtered=length(table(motif_filtered$V4))
  motif_filtered=motif_filtered[!(motif_filtered$V5==-22 & motif_filtered$V6==-1),]
  #allow a slip of 1 bp call methylation
  motif_filtered_filtered <- motif_filtered[motif_filtered$V1 == "GATC" & motif_filtered$V5 == -13 & motif_filtered$V6 == -10 | motif_filtered$V1 == "GATC" & motif_filtered$V5 == -14 & motif_filtered$V6 == -11  | motif_filtered$V1 != "GATC", ]
  motif_filtered_filtered <- motif_filtered_filtered[motif_filtered_filtered$V1 == "CCWGG" & motif_filtered_filtered$V5 == -13 & motif_filtered_filtered$V6 == -9 | motif_filtered_filtered$V1 == "CCWGG" & motif_filtered_filtered$V5 == -15 & motif_filtered_filtered$V6 == -11 | motif_filtered_filtered$V1 != "CCWGG", ]
  motif_filtered_filtered <- motif_filtered_filtered[motif_filtered_filtered$V1 == "AACNNNNNNGTGC" & motif_filtered_filtered$V5 == -13 & motif_filtered_filtered$V6 == -1 | motif_filtered_filtered$V1 == "AACNNNNNNGTGC" & motif_filtered_filtered$V5 == -21 & motif_filtered_filtered$V6 == -9 |motif_filtered_filtered$V1 != "AACNNNNNNGTGC", ]
  motif_filtered_filtered <- motif_filtered_filtered[motif_filtered_filtered$V1 == "GCACNNNNNNGTT" & motif_filtered_filtered$V5 == -14 & motif_filtered_filtered$V6 == -2 | motif_filtered_filtered$V1 == "GCACNNNNNNGTT" & motif_filtered_filtered$V5 == -23 & motif_filtered_filtered$V6 == -11| motif_filtered_filtered$V1 != "GCACNNNNNNGTT", ]
  motif_filtered_filtered <- motif_filtered_filtered[motif_filtered_filtered$V1 == "ATGCAT" & motif_filtered_filtered$V5 == -16 & motif_filtered_filtered$V6 == -11 | motif_filtered_filtered$V1 == "ATGCAT" & motif_filtered_filtered$V5 == -13 & motif_filtered_filtered$V6 == -8 & motif_filtered_filtered$V8==1  | motif_filtered_filtered$V1 != "ATGCAT", ]
  #list of pairs of indices that are duplicated
  duplicated_indices=which(duplicated(motif_filtered_filtered$V4))
  result_list <- list()
  for (num in duplicated_indices) {
    result_vector <- c(num, num - 1)
    result_list <- append(result_list, list(result_vector))
  }
  removal_indices <- c()
  # Iterate through the paired indices and compare scores in V8
  for (pair in result_list) {
    index1 <- pair[1]
    index2 <- pair[2]
    score1 <- as.numeric(motif_filtered_filtered[index1, "V8"])
    score2 <- as.numeric(motif_filtered_filtered[index2, "V8"])
    # Identify the index of the row with the highest score
    if (score1 > score2) {
      removal_indices <- c(removal_indices, index2)
    } else {
      removal_indices <- c(removal_indices, index1)
    }
  }
  # Remove the rows with the indices from the removal list
  data_filtered <- motif_filtered_filtered[-removal_indices, ]
  motif_number_filtered=dim(data_filtered)[1]
  percentage_motifs_filtered=motif_number_filtered/methylation_number_filtered*100
  percentage_motifs_GATC_filtered=sum(data_filtered$V1=="GATC")/motif_number_filtered*100
  percentage_motifs_AACNNNNNNGTGC_filtered=(sum(data_filtered$V1=="AACNNNNNNGTGC")+sum(data_filtered$V1=="GCACNNNNNNGTT"))/motif_number_filtered*100
  percentage_motifs_CCWGG_filtered=sum(data_filtered$V1=="CCWGG")/motif_number_filtered*100
  percentage_motifs_ATGCAT_filtered=sum(data_filtered$V1=="ATGCAT")/motif_number_filtered*100
  score_GATC_filtered=mean(unlist(filter(data_filtered, V1 == "GATC")%>%select(V8)))
  score_AACNNNNNNGTGC_filtered=mean(unlist(filter(data_filtered, V1 == "AACNNNNNNGTGC" | V1 == "GCACNNNNNNGTT")%>%select(V8)))
  score_CCWGG_filtered=mean(unlist(filter(data_filtered, V1 == "CCWGG")%>%select(V8)))
  score_ATGCAT_filtered=mean(unlist(filter(data_filtered, V1 == "ATGCAT")%>%select(V8)))
  
  fraction_partial=methylation_number_filtered/methylation_number_full*100
  fraction_GATC=sum(data_filtered$V1=="GATC")/sum(data_full$V1=="GATC")*100
  fraction_AACNNNNNNGTGC=(sum(data_filtered$V1=="AACNNNNNNGTGC")+sum(data_filtered$V1=="GCACNNNNNNGTT"))/(sum(data_full$V1=="AACNNNNNNGTGC")+sum(data_full$V1=="GCACNNNNNNGTT"))*100
  fraction_CCWGG=sum(data_filtered$V1=="CCWGG")/sum(data_full$V1=="CCWGG")*100
  fraction_ATGCAT=sum(data_filtered$V1=="ATGCAT")/sum(data_full$V1=="ATGCAT")*100
  
  results_df_full=t(data.frame(percentage_motifs_full,percentage_motifs_CCWGG_full,score_CCWGG_full,percentage_motifs_AACNNNNNNGTGC_full,
                               score_AACNNNNNNGTGC_full,percentage_motifs_GATC_full,score_GATC_full,percentage_motifs_ATGCAT_full,score_ATGCAT_full))
  results_df_filtered=t(data.frame(percentage_motifs_filtered,percentage_motifs_CCWGG_filtered,score_CCWGG_filtered,percentage_motifs_AACNNNNNNGTGC_filtered,
                                   score_AACNNNNNNGTGC_filtered,percentage_motifs_GATC_filtered,score_GATC_filtered,percentage_motifs_ATGCAT_filtered,score_ATGCAT_filtered))
  summary_results=t(data.frame(fraction_partial,fraction_GATC,fraction_AACNNNNNNGTGC,fraction_CCWGG,fraction_ATGCAT))
  
  
  #Extract from the fasta file the sequences that have not found a motif
  full_sequences <- readDNAStringSet(methylation_full_fasta)
  df <- data.frame(Seq_Name = names(full_sequences), Sequence = as.character(full_sequences), stringsAsFactors = FALSE)
  sequences=data_full$V4
  sequence_without_motif <- df[!df$Seq_Name %in% sequences, ]
  rownames(sequence_without_motif)=NULL
  dna_strings <- DNAStringSet(sequence_without_motif$Sequence)
  names(dna_strings) <- sequence_without_motif$Seq_Name
  writeXStringSet(dna_strings, file = paste0("no_motif_full_",str_remove(str_remove(motif_22_full,"methylation_"),"_seq_22_full.txt"),".fa"))
  return(list(results_df_full,results_df_filtered,summary_results))
}


results_0_0=motifs(motif_22_filtered = "methylation_0-0_seq_22_filtered.txt",motif_22_full = "methylation_0-0_seq_22_full.txt",methylation_full_fasta = "methylation_0-0_seq_22_full.fa")
results_4_4=motifs(motif_22_filtered = "methylations_4-4_seq_22_filtered.txt",motif_22_full = "methylation_4-4_seq_22_full.txt",methylation_full_fasta = "methylation_4-4_seq_22_full.fa")
results_6_6=motifs(motif_22_filtered = "methylations_6-6_seq_22_filtered.txt",motif_22_full = "methylation_6-6_seq_22_full.txt",methylation_full_fasta = "methylation_6-6_seq_22_full.fa")
results_7_7=motifs(motif_22_filtered = "methylations_7-7_seq_22_filtered.txt",motif_22_full = "methylation_7-7_seq_22_full.txt",methylation_full_fasta = "methylation_7-7_seq_22_full.fa")
results_SRR1536433=motifs(motif_22_filtered = "SRR1536433_seq_22_filtered.txt",motif_22_full = "SRR1536433_seq_22_full.txt",methylation_full_fasta = "SRR1536433_seq_22_full.fa")

motifs_BSeq=function(motif_22_filtered,motif_22_full,methylation_full_fasta){
  motif_full=read.csv(file=motif_22_full,header = FALSE, sep = "\t",comment=";")
  motif_filtered=read.csv(file=motif_22_filtered,header = FALSE, sep = "\t",comment=";")
  methylation_number_full=length(table(motif_full$V4))
  methylation_number_filtered=length(table(motif_filtered$V4))
  fraction_partial=methylation_number_filtered/methylation_number_full*100
  motif_full=motif_full[!(motif_full$V5==-22 & motif_full$V6==-1),]
  motif_full_filtered <- motif_full[motif_full$V1 == "CCWGG" & motif_full$V5 == -13 & motif_full$V6 == -9 | motif_full$V1 == "CCWGG" & motif_full$V5 == -15 & motif_full$V6 == -11 | motif_full$V1 != "CCWGG", ]
  motif_filtered=motif_filtered[!(motif_filtered$V5==-22 & motif_filtered$V6==-1),]
  motif_filtered_filtered <- motif_filtered[motif_filtered$V1 == "CCWGG" & motif_filtered$V5 == -13 & motif_filtered$V6 == -9 | motif_filtered$V1 == "CCWGG" & motif_filtered$V5 == -15 & motif_filtered$V6 == -11 | motif_filtered$V1 != "CCWGG", ]
  motif_number_filtered=dim(motif_filtered_filtered)[1]
  motif_number_full=dim(motif_full_filtered)[1]
  
  
  full_sequences <- readDNAStringSet(methylation_full_fasta)
  df <- data.frame(Seq_Name = names(full_sequences), Sequence = as.character(full_sequences), stringsAsFactors = FALSE)
  sequences=motif_full_filtered$V4
  sequence_without_motif <- df[!df$Seq_Name %in% sequences, ]
  rownames(sequence_without_motif)=NULL
  dna_strings <- DNAStringSet(sequence_without_motif$Sequence)
  names(dna_strings) <- sequence_without_motif$Seq_Name
  writeXStringSet(dna_strings, file = paste0("no_motif_full_",str_remove(motif_22_full,"_full.txt"),".fa"))
  
  
  percentage_CCWGG_filtered=motif_number_filtered/methylation_number_filtered*100
  score_CCWGG_filtered=mean(unlist(filter(motif_filtered_filtered, V1 == "CCWGG")%>%select(V8)))
  percentage_CCWGG_full=motif_number_full/methylation_number_full*100
  score_CCWGG_full=mean(unlist(filter(motif_full_filtered, V1 == "CCWGG")%>%select(V8)))
  fraction_partial=methylation_number_filtered/methylation_number_full*100
  fraction_CCWGG=sum(motif_filtered_filtered$V1=="CCWGG")/sum(motif_full_filtered$V1=="CCWGG")*100
  results=t(data.frame(fraction_partial,percentage_CCWGG_full,score_CCWGG_full,percentage_CCWGG_filtered,score_CCWGG_filtered))
  return(results)
}


results_BS_1=motifs_BSeq(motif_22_filtered = "BSeq_1_filtered.txt",motif_22_full = "BSeq_1_full.txt",methylation_full_fasta = "BSeq_1_full.fa")
results_BS_2=motifs_BSeq(motif_22_filtered = "BSeq_2_filtered.txt",motif_22_full = "BSeq_2_full.txt",methylation_full_fasta = "BSeq_2_full.fa")



