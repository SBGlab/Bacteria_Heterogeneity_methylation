#
motif_conservation=function(motif_5_filtered,motif_5_full,motif_4_filtered,motif_4_full,motif_14_filtered,motif_14_full,motif_6_filtered,motif_6_full){
  #GATC
  motif_4_full=read.csv(file=motif_4_full,header = FALSE, sep = "\t",comment=";")
  motif_4_full_score=motif_4_full[motif_4_full$V8>0,]
  motif_4_full_match=motif_4_full_score[motif_4_full_score$V6==-1 & motif_4_full_score$V5==-4,]
  mean_motif_4_full=mean(motif_4_full_match$V8,)
  sd(motif_4_full_match$V8,)
  motif_4_filtered=read.csv(file=motif_4_filtered,header = FALSE, sep = "\t",comment=";")
  motif_4_filtered_score=motif_4_filtered[motif_4_filtered$V8>0,]
  motif_4_filtered_match=motif_4_filtered_score[motif_4_filtered_score$V6==-1 & motif_4_filtered_score$V5==-4,]
  mean_motif_4_filtered=mean(motif_4_filtered_match$V8,)
  sd(motif_4_filtered_match$V8,)
  percentaje_GATC=nrow(motif_4_filtered_match)/nrow(motif_4_full_match)
  #GCACNNNNNNGTT/AACNNNNNNGTGC
  motif_14_full=read.csv(file=motif_14_full,header = FALSE, sep = "\t",comment=";")
  motif_14_full_score=motif_14_full[motif_14_full$V8>0,]
  motif_14_full_match=motif_14_full_score[motif_14_full_score$V6==-1 & motif_14_full_score$V5==-13 | motif_14_full_score$V6==-2 & motif_14_full_score$V5==-14,]
  mean_motif_14_full=mean(motif_14_full_match$V8,)
  sd(motif_14_full_match$V8,)
  motif_14_filtered=read.csv(file=motif_14_filtered,header = FALSE, sep = "\t",comment=";")
  motif_14_filtered_score=motif_14_filtered[motif_14_filtered$V8>0,]
  motif_14_filtered_match=motif_14_filtered_score[motif_14_filtered_score$V6==-1 & motif_14_filtered_score$V5==-13 | motif_14_filtered_score$V6==-2 & motif_14_filtered_score$V5==-14,]
  mean_motif_14_filtered=mean(motif_14_filtered_match$V8,)
  sd(motif_14_filtered_match$V8,)
  percentage_AACNNNNNNGTGC=nrow(motif_14_filtered_match)/nrow(motif_14_full_match)
  #CCWGG
  motif_5_full=read.csv(file=motif_5_full,header = FALSE, sep = "\t",comment=";")
  motif_5_full_score=motif_5_full[motif_5_full$V8>0,]
  motif_5_full_match=motif_5_full_score[motif_5_full_score$V6==-1 & motif_5_full_score$V5==-5,]
  mean_motif_5_full=mean(motif_5_full_match$V8,)
  sd(motif_5_full_match$V8,)
  motif_5_filtered=read.csv(file=motif_5_filtered,header = FALSE, sep = "\t",comment=";")
  motif_5_filtered_score=motif_5_filtered[motif_5_filtered$V8>0,]
  motif_5_filtered_match=motif_5_filtered_score[motif_5_filtered_score$V6==-1 & motif_5_filtered_score$V5==-5,]
  mean_motif_5_filtered=mean(motif_5_filtered_match$V8,)
  sd(motif_5_filtered_match$V8,)
  motif_6_filtered=read.csv(file=motif_6_filtered,header = FALSE, sep = "\t",comment=";")
  motif_6_filtered_score=motif_6_filtered[motif_6_filtered$V8>0,]
  motif_6_filtered_match=motif_6_filtered_score[motif_6_filtered_score$V6==-1 & motif_6_filtered_score$V5==-6,]
  mean_motif_6_filtered=mean(motif_6_filtered_match$V8,)
  motif_6_full=read.csv(file=motif_6_full,header = FALSE, sep = "\t",comment=";")
  motif_6_full_score=motif_6_full[motif_6_full$V8>0,]
  motif_6_full_match=motif_6_full_score[motif_6_full_score$V6==-1 & motif_6_full_score$V5==-6,]
  mean_motif_6_full=mean(motif_6_full_match$V8,)  
  percentaje_ATGCAT=nrow(motif_6_filtered_match)/nrow(motif_6_full_match)
  percentaje_CCWGG=nrow(motif_5_filtered_match)/nrow(motif_5_full_match)
  proportion_GATC_AACNNNNNNGTGC_full=nrow(motif_4_full_match)/nrow(motif_14_full_match)
  proportion_GATC_AACNNNNNNGTGC_filtered=nrow(motif_4_filtered_match)/nrow(motif_14_filtered_match)
  result=data.frame(mean_motif_6_filtered,mean_motif_6_full,mean_motif_5_filtered,mean_motif_5_full,mean_motif_14_filtered,mean_motif_14_full,mean_motif_4_filtered,mean_motif_4_full,percentaje_ATGCAT,percentaje_GATC,percentage_AACNNNNNNGTGC,percentaje_CCWGG,proportion_GATC_AACNNNNNNGTGC_full,proportion_GATC_AACNNNNNNGTGC_filtered)
  result[,9:12]=result[,9:12]*100
  return (t(result))
}

motif_conservation_BSeq=function(motif_5_filtered,motif_5_full){
  motif_5_full=read.csv(file=motif_5_full,header = FALSE, sep = "\t",comment=";")
  motif_5_full_score=motif_5_full[motif_5_full$V8>0,]
  motif_5_full_match=motif_5_full_score[motif_5_full_score$V6==-1 & motif_5_full_score$V5==-5,]
  mean_motif_5_full=mean(motif_5_full_match$V8,)
  sd(motif_5_full_match$V8,)
  motif_5_filtered=read.csv(file=motif_5_filtered,header = FALSE, sep = "\t",comment=";")
  motif_5_filtered_score=motif_5_filtered[motif_5_filtered$V8>0,]
  motif_5_filtered_match=motif_5_filtered_score[motif_5_filtered_score$V6==-1 & motif_5_filtered_score$V5==-5,]
  mean_motif_5_filtered=mean(motif_5_filtered_match$V8,)
  sd(motif_5_filtered_match$V8,)
  percentaje_CCWGG=nrow(motif_5_filtered_match)/nrow(motif_5_full_match)
  result=data.frame(mean_motif_5_filtered,mean_motif_5_full)
  return (t(result))
}
  

SRR1536433=motif_conservation(motif_6_full="SRR1536433_seq_6_full.txt",motif_6_filtered="SRR1536433_seq_6_filtered.txt",motif_5_full="SRR1536433_seq_5_full.txt",motif_5_filtered="SRR1536433_seq_5_filtered.txt",motif_4_filtered="SRR1536433_seq_4_filtered.txt",motif_4_full="SRR1536433_seq_4_full.txt",motif_14_filtered="SRR1536433_seq_14_filtered.txt",motif_14_full="SRR1536433_seq_14_full.txt")
methylation_0_0=motif_conservation(motif_6_full="methylation_0-0_seq_6_full.txt",motif_6_filtered="methylations_0-0_seq_6_filtered.txt",motif_5_full="methylation_0-0_seq_5_full.txt",motif_5_filtered="methylation_0-0_seq_5_filtered.txt",motif_4_filtered="methylation_0-0_seq_4_filtered.txt",motif_4_full="methylation_0-0_seq_4_full.txt",motif_14_filtered="methylation_0-0_seq_14_filtered.txt",motif_14_full="methylation_0-0_seq_14_full.txt")
methylation_7_7=motif_conservation(motif_6_full="methylation_7-7_seq_6_full.txt",motif_6_filtered = "methylation_7-7_seq_6_filtered.txt",motif_5_full="methylation_7-7_seq_5_full.txt",motif_5_filtered="methylation_7-7_seq_5_filtered.txt",motif_4_filtered="methylation_7-7_seq_4_filtered.txt",motif_4_full="methylation_7-7_seq_4_full.txt",motif_14_filtered="methylation_7-7_seq_14_filtered.txt",motif_14_full="methylation_7-7_seq_14_full.txt")
Bseq=motif_conservation_BSeq(motif_5_filtered="SRP329794_BSeq_seq_5_filtered.txt",motif_5_full="SRP329794_BSeq_seq_5_full.txt")

filtered_5 <- c(SRR1536433[1,], methylation_0_0[1,],methylation_7_7[1,])
full_5 <- c(SRR1536433[2,], methylation_0_0[2,],methylation_7_7[2,])
t_test_5 <- t.test(filtered_5, full_5, paired=TRUE)
filtered_4 <- c(SRR1536433[5,], methylation_0_0[5,],methylation_7_7[5,])
full_4 =c(SRR1536433[6,], methylation_0_0[6,],methylation_7_7[6,])
t_test_4 <- t.test(filtered_4, full_4, paired=TRUE)
filtered_14 <- c(SRR1536433[3,], methylation_0_0[3,],methylation_7_7[3,])
full_14 =c(SRR1536433[4,], methylation_0_0[4,],methylation_7_7[4,])
t_test_14 <- t.test(filtered_14, full_14, paired=TRUE)


