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
  locus_names=names(which(table(pb_annotation$geneId)>=2))
  symbol_names=names(which(table(pb_annotation$transcriptId)>=2))
  matches=data.frame()
  for (i in locus_names){
    matches=rbind(matches,pb_annotation[grep(i,pb_annotation$geneId),])
  }
  write.table(locus_names,file=paste0(str_remove(peaks,"_out.gff"),"_genelist.txt"),row.names=FALSE, col.names=FALSE, sep="",quote=FALSE)
  write.table(matches,file=paste0(str_remove(peaks,"_out.gff"),"_table.txt"),row.names=FALSE, sep="\t",quote=FALSE)
  return (list(locus_names,symbol_names))
}

methylation_0_0_genes=genes_found_locustag_table("methylation_0-0_out.gff")
SRR1536433_genes=genes_found_locustag_table("SRR1536433_out.gff")
SRP329794_BSeq_genes=genes_found_locustag_table("SRP329794_BSeq_out.gff")
methylation_7_7_genes=genes_found_locustag_table("methylation_7-7_out.gff")

functional_enrichment=function(gene_symbol,background=metabolic_genes){
  ego <- enrichGO(gene_symbol,keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = background)
  return(ego)
}

genes_0_0_locus=as.vector(methylation_0_0_genes)[[1]]
genes_7_7_locus=as.vector(methylation_7_7_genes)[[1]]
genes_SRR1536433_locus=as.vector(SRR1536433_genes)[[1]]
genes_Bseq_locus=as.vector(SRP329794_BSeq_genes)[[1]]
genes_0_0_symbol=as.vector(methylation_0_0_genes)[[2]]
genes_7_7_symbol=as.vector(methylation_7_7_genes)[[2]]
genes_SRR1536433_symbol=as.vector(SRR1536433_genes)[[2]]
genes_Bseq_symbol=as.vector(SRP329794_BSeq_genes)[[2]]

functional_0_0=functional_enrichment(genes_0_0_symbol)
barplot(functional_0_0, showCategory=10,title="functional_0_0_gene")
functional_0_0_summary=data.frame(functional_0_0@result)
functional_7_7=functional_enrichment(genes_7_7_symbol)
barplot(functional_7_7, showCategory=10,title="functional_7_7_gene")
functional_7_7_summary=data.frame(functional_7_7@result)
functional_SRR1536433=functional_enrichment(genes_SRR1536433_symbol)
barplot(functional_SRR1536433, showCategory=10,title="functional_SRR1536433_gene")
functional_SRR1536433_summary=data.frame(functional_SRR1536433@result)
functional_Bseq=functional_enrichment(genes_Bseq_symbol)
barplot(functional_Bseq, showCategory=10,title="functional_Bseq_gene")
functional_Bseq_summary=data.frame(functional_Bseq@result)

metabolic_genes_locus=na.exclude(AnnotationDbi::select(txdb, keys=metabolic_genes, keytype="TXNAME",columns="GENEID",multiVals="first")[,2])

KEGG_0_0 <- enrichKEGG(genes_0_0_locus,keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_0_0,showCategory=10,title="KEGG 0_0 gene")
KEGG_7_7 <- enrichKEGG(genes_7_7_locus,keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_7_7,showCategory=10,title="KEGG 7_7 gene")
KEGG_Bseq <- enrichKEGG(genes_Bseq_locus,keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_Bseq,showCategory=10,title="KEGG Bseq gene")
KEGG_SRR1536433 <- enrichKEGG(genes_SRR1536433_locus,keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_SRR1536433,showCategory=10,title="KEGG SRR1536433 gene")


#Genes related with heterogenity in biomass growth and regulated by methylation
expressed_biomass=read.csv("expressed_biomass.txt",header = FALSE, sep = "\n")
expressed_biomass_0_0=expressed_biomass[21:26,]
expressed_biomass_SRR1536433=expressed_biomass[2:5,]
expressed_biomass_7_7=expressed_biomass[7:19,]
expressed_biomass_Bseq=expressed_biomass[29:33,]

locus_to_symbol=function(locus_list){
  return(na.exclude(AnnotationDbi::select(txdb, keys=locus_list, keytype="GENEID",columns="TXNAME",multiVals="first")[,2]))
}

expressed_biomass_0_0_symbol=locus_to_symbol(expressed_biomass_0_0)
expressed_biomass_SRR1536433_symbol=locus_to_symbol(expressed_biomass_SRR1536433)
expressed_biomass_7_7_symbol=locus_to_symbol(expressed_biomass_7_7)
expressed_biomass_Bseq_symbol=locus_to_symbol(expressed_biomass_Bseq)

Bseq_biomass_FE=functional_enrichment(expressed_biomass_Bseq_symbol)
barplot(Bseq_biomass_FE, showCategory=10,title="BSeq biomass gene")
expressed_0_0_biomass_FE=functional_enrichment(expressed_biomass_0_0_symbol)
barplot(expressed_0_0_biomass_FE, showCategory=10,title="expressed_0_0 biomass gene")
expressed_7_7_biomass_FE=functional_enrichment(expressed_biomass_7_7_symbol)
barplot(expressed_7_7_biomass_FE, showCategory=10,title="expressed_7_7 biomass gene")
SRR1536433_biomass_FE=functional_enrichment(expressed_biomass_SRR1536433_symbol)
barplot(SRR1536433_biomass_FE, showCategory=10,title="SRR1536433 biomass gene")

kegg_enrichment=function(locus_list){
  return(enrichKEGG(locus_list,keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus))
}

KEGG_0_0_biomass <- kegg_enrichment(expressed_biomass_0_0)
barplot(KEGG_0_0_biomass,showCategory=10,title="KEGG 0_0 biomass gene")
KEGG_7_7_biomass <- kegg_enrichment(expressed_biomass_7_7)
barplot(KEGG_7_7_biomass,showCategory=10,title="KEGG 7_7  biomass gene")
KEGG_Bseq_biomass <- kegg_enrichment(expressed_biomass_Bseq)
barplot(KEGG_Bseq_biomass,showCategory=10,title="KEGG Bseq  biomass gene")
KEGG_SRR1536433_biomass <- kegg_enrichment(expressed_biomass_SRR1536433)
barplot(KEGG_SRR1536433_biomass,showCategory=10,title="KEGG SRR1536433  biomass gene")


#Genes related with heterogenity NOT biomass growth and regulated by methylation
expressed_not_biomass=read.csv("expressed_not_biomass.txt",header = FALSE, sep = "\n")
expressed_not_biomass_0_0=expressed_not_biomass[27:35,]
expressed_not_biomass_SRR1536433=expressed_not_biomass[2:8,]
expressed_not_biomass_7_7=expressed_not_biomass[10:25,]
expressed_not_biomass_Bseq=expressed_not_biomass[37:42,]

expressed_not_biomass_0_0_symbol=locus_to_symbol(expressed_not_biomass_0_0)
expressed_not_biomass_SRR1536433_symbol=locus_to_symbol(expressed_not_biomass_SRR1536433)
expressed_not_biomass_7_7_symbol=locus_to_symbol(expressed_not_biomass_7_7)
expressed_not_biomass_Bseq_symbol=locus_to_symbol(expressed_not_biomass_Bseq)

Bseq_not_biomass_FE=functional_enrichment(expressed_not_biomass_Bseq_symbol)
barplot(Bseq_not_biomass_FE, showCategory=10,title="BSeq not biomass gene")
expressed_0_0_not_biomass_FE=functional_enrichment(expressed_not_biomass_0_0_symbol)
barplot(expressed_0_0_not_biomass_FE, showCategory=10,title="expressed_0_0 not biomass gene")
expressed_7_7_not_biomass_FE=functional_enrichment(expressed_not_biomass_7_7_symbol)
barplot(expressed_7_7_not_biomass_FE, showCategory=10,title="expressed_7_7 not biomass gene")
SRR1536433_not_biomass_FE=functional_enrichment(expressed_not_biomass_SRR1536433_symbol)
barplot(SRR1536433_not_biomass_FE, showCategory=10,title="SRR1536433 not biomass gene")

KEGG_0_0_not_biomass <- kegg_enrichment(expressed_not_biomass_0_0)
barplot(KEGG_0_0_not_biomass,showCategory=10,title="KEGG 0_0 not biomass gene")
KEGG_7_7_not_biomass <- kegg_enrichment(expressed_not_biomass_7_7)
barplot(KEGG_7_7_not_biomass,showCategory=10,title="KEGG 7_7 not biomass gene")
KEGG_Bseq_not_biomass <- kegg_enrichment(expressed_not_biomass_Bseq)
barplot(KEGG_Bseq_not_biomass,showCategory=10,title="KEGG Bseq not biomass gene")
KEGG_SRR1536433_not_biomass <- kegg_enrichment(expressed_not_biomass_SRR1536433)
barplot(KEGG_SRR1536433_not_biomass,showCategory=10,title="KEGG SRR1536433 not biomass gene")
