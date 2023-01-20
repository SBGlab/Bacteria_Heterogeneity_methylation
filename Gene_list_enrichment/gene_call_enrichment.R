library(GenomicFeatures)
library(org.EcK12.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(dplyr)
library(stringr)
txdb <- makeTxDbFromGFF("GCF_000005845.2_ASM584v2_genomic.gff")
metabolic_genes=read.csv("conv_38DC2E424FE51674044440701.txt",sep="\t",colClasses = c("NULL","character","NULL","NULL"))
metabolic_genes=as.vector(metabolic_genes)[[1]]

genes_found_locustag_table=function(peaks){
  peakAnnoList2 <- annotatePeak(peaks,TxDb=txdb,tssRegion=c(-1000,0),verbose=TRUE,overlap="all")
  pb_annotation2 <- data.frame(peakAnnoList2@anno)
  peakAnnoList3 <- annotatePeak(peaks,TxDb=txdb,tssRegion=c(-1000,0),verbose=FALSE)
  pb_annotation3 <- data.frame(peakAnnoList3@anno)
  pb_annotation=rbind(pb_annotation2,pb_annotation3)
  pb_annotation=pb_annotation %>% distinct()
  pb_annotation=pb_annotation[order(pb_annotation$start),]
  pb_annotation$geneStrand[pb_annotation$geneStrand == 1] <- '+'
  pb_annotation$geneStrand[pb_annotation$geneStrand == 2] <- '-'
  pb_annotation=pb_annotation[pb_annotation$STRAND == pb_annotation$geneStrand,]
  gene_names=names(which(table(pb_annotation$geneId)>=2))
  matches=data.frame()
  for (i in gene_names){
    matches=rbind(matches,pb_annotation[grep(i,pb_annotation$geneId),])
  }
  write.table(gene_names,file=paste0(str_remove(peaks,"_out.gff"),"_genelist.txt"),row.names=FALSE, col.names=FALSE, sep="",quote=FALSE)
  write.table(matches,file=paste0(str_remove(peaks,"_out.gff"),"_table.txt"),row.names=FALSE, sep="\t",quote=FALSE)
  return (list(gene_names,matches))
}

methylation_0_0_genes=genes_found_locustag_table("methylation_0-0_out.gff")
SRR1536433_genes=genes_found_locustag_table("SRR1536433_out.gff")
SRP329794_BSeq_genes=genes_found_locustag_table("SRP329794_BSeq_out.gff")
methylation_7_7_genes=genes_found_locustag_table("methylation_7-7_out.gff")


functional_enrichment=function(peaks,metabolic_genes){
  peakAnnoList2 <- annotatePeak(peaks,TxDb=txdb,tssRegion=c(-1000,0),verbose=TRUE,overlap="all")
  pb_annotation2 <- data.frame(peakAnnoList2@anno)
  peakAnnoList3 <- annotatePeak(peaks,TxDb=txdb,tssRegion=c(-1000,0),verbose=FALSE)
  pb_annotation3 <- data.frame(peakAnnoList3@anno)
  pb_annotation=rbind(pb_annotation2,pb_annotation3)
  pb_annotation=pb_annotation %>% distinct()
  pb_annotation=pb_annotation[order(pb_annotation$start),]
  pb_annotation$geneStrand[pb_annotation$geneStrand == 1] <- '+'
  pb_annotation$geneStrand[pb_annotation$geneStrand == 2] <- '-'
  pb_annotation=pb_annotation[pb_annotation$STRAND == pb_annotation$geneStrand,]
  gene_names=names(which(table(pb_annotation$transcriptId)>=2))
  ego <- enrichGO(as.vector(gene_names),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe =metabolic_genes)
  cluster_summary=data.frame(ego@result)
  write.table(cluster_summary,file=paste0(str_remove(peaks,"_out.gff"),"_FE.txt"),row.names = FALSE,sep="\t")
  return (cluster_summary)
}


methylation_0_0_FE=functional_enrichment("methylation_0-0_out.gff",metabolic_genes)
methylation_7_7_FE=functional_enrichment("methylation_7-7_out.gff",metabolic_genes)
SRR1536433_FE=functional_enrichment("SRR1536433_out.gff")
SRP329794_BSeq_FE=functional_enrichment("SRP329794_BSeq_out.gff")

