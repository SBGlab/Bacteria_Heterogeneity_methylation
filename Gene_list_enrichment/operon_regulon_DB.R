library(GenomicFeatures)
library(dplyr)
library(stringr)
library(tidyverse)
library(gridExtra)
library(clusterProfiler)
library(org.EcK12.eg.db)
setwd("E:/Universidad/MBB/TFM/gene_call_enrichment/operon")

operon=read.csv("operon.txt",header = FALSE, sep = "\n")
mylist_split <- lapply(operon[39:nrow(operon),], function(x) strsplit(x, "\t"))
x=t(as.data.frame(mylist_split))
x=x[,2:7]
x[,6][x[,6] == "forward"] <- '+'
x[,6][x[,6] == "reverse"] <- '-'
x=as.data.frame(x)
#add an extra column with the chromosome name to convert to bed format
x <- data.frame(chromosome = 'NC_000913', x)
colnames(x)=c("chromosome","operon name", "first genepos left","last geneposright","regulation pos left","regulation pos right","operon strand")
rownames(x)=NULL

#These genes are interrumped by a mobile element-->eliminate
#NC_000913	2169693	2171727	-	gatR
#NC_000913	2068635	2070474	+	yoeA
#NC_000913	390251	394418	+	yaiT
#NC_000913	3651291	3653713	+	yhiS
#NC_000913	380844	382872	-	yaiX
#NC_000913	1465392	1474013	+	ydbA
x=x[x[,4]!=0,]

operon_reg <- data.frame()
for (i in 1:nrow(x)) {
  strand <- x[i, 7]
  if (strand == "+") {
    new_row <- as.vector(unlist(x[i, c(1, 2, 5,3, 7)]))
    trail_regulatory=as.vector(unlist(x[i, c(1, 2, 4,6, 7)]))
  } 
  else if (strand == "-") {
    new_row <- as.vector(unlist(x[i, c(1, 2, 4, 6, 7), drop = FALSE]))
    trail_regulatory=as.vector(unlist(x[i, c(1, 2, 5,3, 7)]))
  }
  if (trail_regulatory[3]!=trail_regulatory[4]){
    operon_reg <- rbind(operon_reg, new_row,trail_regulatory)
  }
  else{
    operon_reg <- rbind(operon_reg, new_row)
  }
  colnames(operon_reg)=NULL
}
operon_reg=operon_reg[, c(1,3,4,5,2)]
operon_reg=operon_reg[!duplicated(operon_reg), ]
##Eliminate regulatory regions lenght 0
operon_reg=operon_reg[operon_reg[,2] != operon_reg[,3],]
#mean of size of regulatory region
mean(as.integer(operon_reg[,3])-as.integer(operon_reg[,2]))

operon_cds=data.frame()
for (i in 1:nrow(x)) {
  strand <- x[i, 7]
  if (strand == "+") {
    new_row <- as.vector(unlist(x[i, c(1, 2, 3, 4, 7)]))
  }
  else if (strand == "-") {
    new_row <- as.vector(unlist(x[i, c(1, 2,3,4, 7), drop = FALSE]))
  }
  operon_cds <- rbind(operon_cds, new_row)
  colnames(operon_cds)=NULL
}
operon_cds=operon_cds[, c(1,3,4,5,2)]
write.table(operon_cds,file="operon_cds.bed",col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")
write.table(operon_reg,file="operon_reg.bed",col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")

###bedtools intersect (server)

operon_gene=read.csv("OperonSet.txt",sep="\t",comment = "#",header=TRUE)[,c(1,6)]
txdb <- makeTxDbFromGFF("GCF_000005845.2_ASM584v2_genomic.gff")

#Order of significance of operons and associated genes
significance=function(file_cds,file_reg,threshold_cds,threshold_reg){
  operon_cds=read.csv(file = file_cds,header = FALSE, sep = "\t")[,10]
  operon_order_cds=data.frame(sort(table(operon_cds),decreasing = TRUE))
  operon_meth_cds=operon_order_cds[operon_order_cds[,2]>=threshold_cds & !operon_order_cds[,1]=="micA" & !operon_order_cds[,1]=="dnaX" & !operon_order_cds[,1]=="yneO" & !operon_order_cds[,1]=="copA",]
  operon_reg=read.csv(file = file_reg,header = FALSE, sep = "\t")[,10]
  operon_order_reg=data.frame(sort(table(operon_reg),decreasing = TRUE))
  operon_meth_reg=operon_order_reg[operon_order_reg[,2]>=threshold_reg & !(operon_order_reg[,1]=="dnaX") & !(operon_order_reg[,1]=="micA") & !(operon_order_reg[,1]=="proK"),]
  colnames(operon_meth_reg)=c("operon","frequency")
  colnames(operon_meth_cds)=c("operon","frequency")
  joined_df <- full_join(operon_meth_reg, operon_meth_cds, by = "operon")
  joined_df$frequency.x <- ifelse(!is.na(joined_df$frequency.x), joined_df$frequency.x, 0)
  joined_df$frequency.y<- ifelse(!is.na(joined_df$frequency.y), joined_df$frequency.y, 0)
  colnames(joined_df)=c("operon","Frec reg","Frec cds")
  for (name in 1:nrow(joined_df)){
    gene_symbol_vector=operon_gene[(operon_gene[,1] %in% joined_df[name,1]),2]
    genes_symbol=unlist(strsplit(gene_symbol_vector,","))
    converted_names=na.exclude(AnnotationDbi::select(txdb, keys=genes_symbol, keytype="TXNAME",columns="GENEID",multiVals="first")[,2])
    result <- paste(converted_names, collapse = ", ")
    joined_df[name,4]=result
    joined_df[name,5]=gene_symbol_vector
    colnames(joined_df)[c(4,5)]=c("Locus tag","Gene symbol")
  }
  gene_names=unlist(strsplit(joined_df$`Locus tag`,", "))
  write.table(gene_names,file=paste0(str_remove(str_remove(file_cds,"overlaps_"),"_cds.txt"),"_genelist.txt"),row.names=FALSE, col.names=FALSE, sep="",quote=FALSE)
  return(joined_df)
}

table_Bseq_1=significance(file_cds = "operon_overlaps_BS_1_cds.txt",file_reg = "operon_overlaps_BS_1_reg.txt",threshold_cds = 4,threshold_reg = 1)
table_Bseq_2=significance(file_cds = "operon_overlaps_BS_2_cds.txt",file_reg = "operon_overlaps_BS_2_reg.txt",threshold_cds = 4,threshold_reg = 1)
table_0_0=significance(file_cds = "operon_overlaps_0_0_cds.txt",file_reg = "operon_overlaps_0_0_reg.txt",threshold_cds = 4,threshold_reg = 1)
table_7_7=significance(file_cds = "operon_overlaps_7_7_cds.txt",file_reg = "operon_overlaps_7_7_reg.txt",threshold_cds = 4,threshold_reg = 1)
table_6_6=significance(file_cds = "operon_overlaps_6_6_cds.txt",file_reg = "operon_overlaps_6_6_reg.txt",threshold_cds = 4,threshold_reg = 1)
table_4_4=significance(file_cds = "operon_overlaps_4_4_cds.txt",file_reg = "operon_overlaps_4_4_reg.txt",threshold_cds = 4,threshold_reg = 1)
table_SRR1536433=significance(file_cds = "operon_overlaps_SRR1536433_cds.txt",file_reg = "operon_overlaps_SRR1536433_reg.txt",threshold_cds = 4,threshold_reg = 1)
genes_SRR1536433_locus_all=unlist(strsplit(table_SRR1536433$`Locus tag`,", "))
genes_0_0_locus_all= unlist(strsplit(table_0_0$`Locus tag`,", "))
genes_6_6_locus_all= unlist(strsplit(table_6_6$`Locus tag`,", "))
genes_Bseq_1_locus_all=unlist(strsplit(table_Bseq_1$`Locus tag`,", "))
genes_Bseq_2_locus_all=unlist(strsplit(table_Bseq_2$`Locus tag`,", "))
genes_7_7_locus_all=unlist(strsplit(table_7_7$`Locus tag`,", "))
genes_4_4_locus_all=unlist(strsplit(table_4_4$`Locus tag`,", "))
genes_0_0_symbol_all=unlist(strsplit(table_0_0$`Gene symbol`,","))
genes_7_7_symbol_all=unlist(strsplit(table_7_7$`Gene symbol`,","))
genes_4_4_symbol_all=unlist(strsplit(table_4_4$`Gene symbol`,","))
genes_6_6_symbol_all=unlist(strsplit(table_6_6$`Gene symbol`,","))
genes_SRR1536433_symbol_all=unlist(strsplit(table_SRR1536433$`Gene symbol`,","))
genes_Bseq_1_symbol_all=unlist(strsplit(table_Bseq_1$`Gene symbol`,","))
genes_Bseq_2_symbol_all=unlist(strsplit(table_Bseq_2$`Gene symbol`,","))


sample_table_7_7=table_7_7[sort(sample(nrow(table_7_7), 10)),]
write.table(sample_table_7_7,file="sample_table_7_7.txt",sep="\t",quote = FALSE)
#Gene intersections all

all_intersections <- list()
for (i in 2:7) {
  combos <- combn(list(genes_SRR1536433_locus_all, genes_0_0_locus_all,genes_Bseq_1_locus_all,genes_Bseq_2_locus_all,genes_7_7_locus_all,genes_4_4_locus_all,genes_6_6_locus_all), i, simplify = FALSE)
  for (j in combos) {
    int <- Reduce(intersect, j)
    #This is the same as doing (for i=4) intersect(unlist(combos[[1]][1]),unlist(combos[[1]][2])) only for two vectors possible
    names_j <- names(list(genes_SRR1536433=genes_SRR1536433_locus_all, genes_0_0=genes_0_0_locus_all, genes_Bseq_1=genes_Bseq_1_locus_all, genes_Bseq_2=genes_Bseq_2_locus_all,genes_7_7=genes_7_7_locus_all,genes_4_4=genes_4_4_locus_all,genes_6_6=genes_6_6_locus_all))[match(j, list(genes_SRR1536433=genes_SRR1536433_locus_all, genes_0_0=genes_0_0_locus_all, genes_Bseq_1=genes_Bseq_1_locus_all, genes_Bseq_2=genes_Bseq_2_locus_all, genes_7_7=genes_7_7_locus_all,genes_4_4=genes_4_4_locus_all,genes_6_6=genes_6_6_locus_all))]
    name <- paste(names_j, collapse = " vs ")
    all_intersections[[name]] <- int
  }
}

all_intersections


metabolic_genes=read.csv("conv_38DC2E424FE51674044440701.txt",sep="\t",colClasses = c("character","character","NULL","NULL"))
metabolic_genes_locus=as.vector(metabolic_genes)[[1]]
metabolic_genes_symbol=as.vector(metabolic_genes)[[2]]


#from the table keep only the rows that have metabolic genes 

gene_list_metabolic=function(table){
  process_locus <- function(row1) {
    row1[4] <- paste(intersect(strsplit(row1[4], ", ")[[1]], metabolic_genes_locus), collapse = ",")
    return(row1)
  }
  process_symbol <- function(row2) {
    row2[5] <- paste(intersect(strsplit(row2[5], ",")[[1]], metabolic_genes_symbol), collapse = ",")
    return(row2)
  }
  table_metabolic <- as.data.frame(t(apply(table, 1, process_locus)))
  table_metabolic <- table_metabolic[table_metabolic[,4] != "",]
  table_metabolic <- as.data.frame(t(apply(table_metabolic, 1, process_symbol)))
  table_metabolic <- table_metabolic[table_metabolic[,5] != "",]
  gene_names=unlist(strsplit(table_metabolic$`Locus tag`,", "))
  write.table(gene_names,file=paste0("operon_",str_remove(deparse(substitute(table)),"table_"),"_genelist_metabolic.txt"),row.names=FALSE, col.names=FALSE, sep="",quote=FALSE)
  return(table_metabolic)
}



table_Bseq_1_metabolic=gene_list_metabolic(table_Bseq_1)
table_Bseq_2_metabolic=gene_list_metabolic(table_Bseq_2)
table_0_0_metabolic=gene_list_metabolic(table_0_0)
table_7_7_metabolic=gene_list_metabolic(table_7_7)
table_SRR1536433_metabolic=gene_list_metabolic(table_SRR1536433)
table_4_4_metabolic=gene_list_metabolic(table_4_4)
table_6_6_metabolic=gene_list_metabolic(table_6_6)


#Fisher exact test to evaluate enrichment on metabolic genes
fisher_test_metabolic=function(table_metabolic,table_all){
  k=dim(table_metabolic)[1]
  n=dim(table_all)[1]
  G=4256
  m=length(metabolic_genes_locus)
  contingency_table <- matrix(c(k, n-k, m-k, G-m-n+k), nrow = 2, byrow = TRUE,dimnames = list(c("Metabolic", "Non-metabolic"),c("My sample", "Out of my sample")))
  return(fisher.test(contingency_table)$p.value)
}

fisher_BSeq_1=fisher_test_metabolic(table_Bseq_1_metabolic,table_Bseq_1)
fisher_BSeq_2=fisher_test_metabolic(table_Bseq_2_metabolic,table_Bseq_2)
fisher_0_0=fisher_test_metabolic(table_0_0_metabolic,table_0_0)
fisher_4_4=fisher_test_metabolic(table_7_7_metabolic,table_7_7)
fisher_6_6=fisher_test_metabolic(table_SRR1536433_metabolic,table_SRR1536433)
fisher_7_7=fisher_test_metabolic(table_4_4_metabolic,table_4_4)
fisher_SRR1536433=fisher_test_metabolic(table_6_6_metabolic,table_6_6)

pvalues_fisher_test=format(t(data.frame(fisher_BSeq_1,fisher_BSeq_2,fisher_0_0,fisher_4_4,fisher_6_6,fisher_7_7,fisher_SRR1536433)),digits=3,scientific=FALSE)
colnames(pvalues_fisher_test)=c("pvalue")

genes_SRR1536433_locus_metabolic=unlist(strsplit(table_SRR1536433_metabolic$`Locus tag`,", "))
genes_0_0_locus_metabolic= unlist(strsplit(table_0_0_metabolic$`Locus tag`,", "))
genes_Bseq_1_locus_metabolic=unlist(strsplit(table_Bseq_1_metabolic$`Locus tag`,", "))
genes_Bseq_2_locus_metabolic=unlist(strsplit(table_Bseq_2_metabolic$`Locus tag`,", "))
genes_4_4_locus_metabolic= unlist(strsplit(table_4_4_metabolic$`Locus tag`,", "))
genes_7_7_locus_metabolic=unlist(strsplit(table_7_7_metabolic$`Locus tag`,", "))
genes_6_6_locus_metabolic=unlist(strsplit(table_6_6_metabolic$`Locus tag`,", "))
genes_0_0_symbol_metabolic=unlist(strsplit(table_0_0_metabolic$`Gene symbol`,","))
genes_7_7_symbol_metabolic=unlist(strsplit(table_7_7_metabolic$`Gene symbol`,","))
genes_SRR1536433_symbol_metabolic=unlist(strsplit(table_SRR1536433_metabolic$`Gene symbol`,","))
genes_Bseq_1_symbol_metabolic=unlist(strsplit(table_Bseq_1_metabolic$`Gene symbol`,","))
genes_Bseq_2_symbol_metabolic=unlist(strsplit(table_Bseq_2_metabolic$`Gene symbol`,","))
genes_4_4_symbol_metabolic=unlist(strsplit(table_4_4_metabolic$`Gene symbol`,","))
genes_6_6_symbol_metabolic=unlist(strsplit(table_6_6_metabolic$`Gene symbol`,","))


#Gene intersections

metabolic_intersections <- list()
for (i in 2:7) {
  combos <- combn(list(genes_SRR1536433_locus_metabolic, genes_0_0_locus_metabolic,genes_Bseq_1_locus_metabolic,genes_Bseq_2_locus_metabolic,genes_7_7_locus_metabolic,genes_4_4_locus_metabolic,genes_6_6_locus_metabolic), i, simplify = FALSE)
  for (j in combos) {
    int <- Reduce(intersect, j)
    names_j <- names(list(genes_SRR1536433=genes_SRR1536433_locus_metabolic, genes_0_0=genes_0_0_locus_metabolic, genes_Bseq_1=genes_Bseq_1_locus_metabolic,genes_Bseq_2=genes_Bseq_2_locus_metabolic, genes_7_7=genes_7_7_locus_metabolic,genes_4_4=genes_4_4_locus_metabolic,genes_6_6=genes_6_6_locus_metabolic))[match(j, list(genes_SRR1536433=genes_SRR1536433_locus_metabolic, genes_0_0=genes_0_0_locus_metabolic, genes_Bseq_1=genes_Bseq_1_locus_metabolic,genes_Bseq_2=genes_Bseq_2_locus_metabolic, genes_7_7=genes_7_7_locus_metabolic,genes_4_4=genes_4_4_locus_metabolic,genes_6_6=genes_6_6_locus_metabolic))]
    name <- paste(names_j, collapse = " vs ")
    metabolic_intersections[[name]] <- int
  }
}
metabolic_intersections

#Functional Analysis
#metabolic_genes_locus=na.exclude(AnnotationDbi::select(txdb, keys=metabolic_genes, keytype="TXNAME",columns="GENEID",multiVals="first")[,2])
functional_enrichment_metabolic=function(gene_symbol,background=metabolic_genes_symbol){
  ego <- enrichGO(as.vector(gene_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.0001,readable = TRUE,universe = background)
  return(ego)
}

functional_enrichment_all_background=function(gene_symbol){
  ego <- enrichGO(as.vector(gene_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.0001,readable = TRUE)
  return(ego)
}

functional_0_0_metabolic=functional_enrichment_metabolic(genes_0_0_symbol_metabolic)
barplot(functional_0_0_metabolic, showCategory=20,title="functional_0_0_operon_metabolic")
functional_0_0_summary=data.frame(functional_0_0_metabolic@result)
functional_0_0_all=functional_enrichment_all_background(genes_0_0_symbol_all)
functional_0_0_all_summary=data.frame(functional_0_0_all@result)
barplot(functional_0_0_all, showCategory=10,title="functional_0_0_operon_all")

functional_7_7_metabolic=functional_enrichment_metabolic(genes_7_7_symbol_metabolic)
barplot(functional_7_7_metabolic, showCategory=30,title="functional_7_7_operon_metabolic")
functional_7_7_summary=data.frame(functional_7_7_metabolic@result)
functional_7_7_all=functional_enrichment_all_background(genes_7_7_symbol_all)
barplot(functional_7_7_all, showCategory=10,title="functional_7_7_operon_all")

functional_4_4_metabolic=functional_enrichment_metabolic(genes_4_4_symbol_metabolic)
barplot(functional_4_4_metabolic, showCategory=30,title="functional_4_4_operon_metabolic")
functional_4_4_summary=data.frame(functional_4_4_metabolic@result)
functional_4_4_all=functional_enrichment_all_background(genes_4_4_symbol_all)
barplot(functional_4_4_all, showCategory=10,title="functional_4_4_operon_all")

functional_6_6_metabolic=functional_enrichment_metabolic(genes_6_6_symbol_metabolic)
barplot(functional_6_6_metabolic, showCategory=30,title="functional_6_6_operon_metabolic")
functional_6_6_summary=data.frame(functional_6_6_metabolic@result)
functional_6_6_all=functional_enrichment_all_background(genes_6_6_symbol_all)
barplot(functional_6_6_all, showCategory=10,title="functional_6_6_operon_all")

functional_SRR1536433_metabolic=functional_enrichment_metabolic(genes_SRR1536433_symbol_metabolic)
barplot(functional_SRR1536433_metabolic, showCategory=10,title="functional_SRR1536433_operon_metabolic")
functional_SRR1536433_summary=data.frame(functional_SRR1536433_metabolic@result)
functional_SRR1536433_all=functional_enrichment_all_background(genes_SRR1536433_symbol_all)
barplot(functional_SRR1536433_all, showCategory=10,title="functional_SRR1536433_operon_all")

functional_Bseq_1_metabolic=functional_enrichment_metabolic(genes_Bseq_1_symbol_metabolic)
barplot(functional_Bseq_1_metabolic, showCategory=10,title="functional_Bseq_1_operon_metabolic")
functional_Bseq_1_summary=data.frame(functional_Bseq_1_metabolic@result)
functional_Bseq_1_all=functional_enrichment_all_background(genes_Bseq_1_symbol_all)
barplot(functional_Bseq_1_all, showCategory=10,title="functional_Bseq_1_operon_all")

functional_BSeq_2_metabolic=functional_enrichment_metabolic(genes_Bseq_2_symbol_metabolic)
barplot(functional_BSeq_2_metabolic, showCategory=10,title="functional_BSeq_2_operon_metabolic")
functional_BSeq_2_summary=data.frame(functional_BSeq_2_metabolic@result)
functional_BSeq_2_all=functional_enrichment_all_background(genes_Bseq_2_symbol_all)
barplot(functional_BSeq_2_all, showCategory=10,title="functional_BSeq_2_operon_all")

KEGG_0_0 <- enrichKEGG(as.vector(genes_0_0_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_0_0,showCategory=10,title="KEGG 0_0 operon")
KEGG_6_6 <- enrichKEGG(as.vector(genes_6_6_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_6_6,showCategory=10,title="KEGG 6_6 operon")
KEGG_4_4 <- enrichKEGG(as.vector(genes_4_4_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_4_4,showCategory=10,title="KEGG 4_4 operon")
KEGG_7_7 <- enrichKEGG(as.vector(genes_7_7_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_7_7,showCategory=10,title="KEGG 7_7 operon")
KEGG_Bseq_1 <- enrichKEGG(as.vector(genes_Bseq_1_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_Bseq_1,showCategory=10,title="KEGG Bseq_1 operon")
KEGG_Bseq_2 <- enrichKEGG(as.vector(genes_Bseq_2_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_Bseq_2,showCategory=10,title="KEGG Bseq_2 operon")
KEGG_SRR1536433 <- enrichKEGG(as.vector(genes_SRR1536433_locus_all),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_SRR1536433,showCategory=10,title="KEGG SRR1536433 operon")








##Related to Alvaro's work

#Genes related with heterogenity in biomass growth and regulated by methylation 
expressed_biomass=read.csv("expressed_biomass.txt",header = FALSE, sep = "\n")
expressed_biomass_0_0=expressed_biomass[63:85,]
expressed_biomass_SRR1536433=expressed_biomass[2:9,]
expressed_biomass_7_7=expressed_biomass[11:61,]
expressed_biomass_Bseq=expressed_biomass[87:112,]
expressed_biomass_4_4=expressed_biomass[114:150,]
expressed_biomass_6_6=expressed_biomass[152:196,]
expressed_biomass_0_0_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_0_0, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_SRR1536433_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_SRR1536433, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_7_7_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_7_7, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_Bseq_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_Bseq, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_4_4_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_4_4, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_6_6_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_6_6, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])

Bseq_biomass_FE=enrichGO(as.vector(expressed_biomass_Bseq_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes_symbol)
barplot(Bseq_biomass_FE, showCategory=10,title="BSeq biomass operon")
expressed_0_0_biomass_FE=enrichGO(as.vector(expressed_biomass_0_0_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes_symbol)
barplot(expressed_0_0_biomass_FE, showCategory=10,title="expressed_0_0 biomass operon")
expressed_7_7_biomass_FE=enrichGO(as.vector(expressed_biomass_7_7_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes_symbol)
barplot(expressed_7_7_biomass_FE, showCategory=10,title="expressed_7_7 biomass operon")
SRR1536433_biomass_FE=enrichGO(as.vector(expressed_biomass_SRR1536433_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes_symbol)
barplot(SRR1536433_biomass_FE, showCategory=10,title="SRR1536433 biomass operon")
expressed_4_4_biomass_FE=enrichGO(as.vector(expressed_biomass_4_4_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes_symbol)
barplot(expressed_4_4_biomass_FE, showCategory=10,title="expressed_4_4 biomass operon")
expressed_6_6_biomass_FE=enrichGO(as.vector(expressed_biomass_6_6_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes_symbol)
barplot(expressed_6_6_biomass_FE, showCategory=10,title="expressed_6_6 biomass operon")

KEGG_0_0_biomass <- enrichKEGG(as.vector(expressed_biomass_0_0),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_0_0_biomass,showCategory=10,title="KEGG 0_0 biomass operon")
KEGG_7_7_biomass <- enrichKEGG(as.vector(expressed_biomass_7_7),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_7_7_biomass,showCategory=10,title="KEGG 7_7  biomass operon")
KEGG_Bseq_biomass <- enrichKEGG(as.vector(expressed_biomass_Bseq),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_Bseq_biomass,showCategory=10,title="KEGG Bseq  biomass operon")
KEGG_SRR1536433_biomass <- enrichKEGG(as.vector(expressed_biomass_SRR1536433),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_SRR1536433_biomass,showCategory=10,title="KEGG SRR1536433  biomass operon")


#Genes related with heterogenity NOT biomass growth and regulated by methylation
expressed_not_biomass=read.csv("expressed_not_biomass.txt",header = FALSE, sep = "\n")
expressed_not_biomass_0_0=expressed_not_biomass[98:126,]
expressed_not_biomass_SRR1536433=expressed_not_biomass[2:18,]
expressed_not_biomass_7_7=expressed_not_biomass[20:96,]
expressed_not_biomass_Bseq=expressed_not_biomass[128:138,]
expressed_not_biomass_0_0_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_not_biomass_0_0, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_not_biomass_SRR1536433_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_not_biomass_SRR1536433, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_not_biomass_7_7_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_not_biomass_7_7, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_not_biomass_Bseq_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_not_biomass_Bseq, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
Bseq_not_biomass_FE=enrichGO(as.vector(expressed_not_biomass_Bseq_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(Bseq_not_biomass_FE, showCategory=10,title="BSeq not biomass operon")
expressed_0_0_not_biomass_FE=enrichGO(as.vector(expressed_not_biomass_0_0_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(expressed_0_0_not_biomass_FE, showCategory=10,title="expressed_0_0 not biomass operon")
expressed_7_7_not_biomass_FE=enrichGO(as.vector(expressed_not_biomass_7_7_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(expressed_7_7_not_biomass_FE, showCategory=10,title="expressed_7_7 not biomass operon")
SRR1536433_not_biomass_FE=enrichGO(as.vector(expressed_not_biomass_SRR1536433_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(SRR1536433_not_biomass_FE, showCategory=10,title="SRR1536433 not biomass operon")

KEGG_0_0_not_biomass <- enrichKEGG(as.vector(expressed_not_biomass_0_0),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_0_0_not_biomass,showCategory=10,title="KEGG 0_0 not biomass operon")
KEGG_7_7_not_biomass <- enrichKEGG(as.vector(expressed_not_biomass_7_7),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_7_7_not_biomass,showCategory=10,title="KEGG 7_7 not biomass operon")
KEGG_Bseq_not_biomass <- enrichKEGG(as.vector(expressed_not_biomass_Bseq),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_Bseq_not_biomass,showCategory=10,title="KEGG Bseq not biomass operon")
KEGG_SRR1536433_not_biomass <- enrichKEGG(as.vector(expressed_not_biomass_SRR1536433),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_SRR1536433_not_biomass,showCategory=10,title="KEGG SRR1536433 not biomass operon")

essential_genes=as.vector((read.csv("essential_genes.txt",header = FALSE, sep = "\n")[1]))[[1]]
essential_genes_symbol=na.exclude(AnnotationDbi::select(txdb, keys=essential_genes, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
genes_0_0_symbol_metabolic[genes_0_0_symbol_metabolic %in% essential_genes_symbol]
genes_4_4_symbol_metabolic[genes_4_4_symbol_metabolic %in% essential_genes_symbol]
genes_6_6_symbol_metabolic[genes_6_6_symbol_metabolic %in% essential_genes_symbol]
genes_7_7_symbol_metabolic[genes_7_7_symbol_metabolic %in% essential_genes_symbol]
genes_SRR1536433_symbol_metabolic[genes_SRR1536433_symbol_metabolic %in% essential_genes_symbol]
genes_Bseq_symbol_metabolic[genes_Bseq_symbol_metabolic %in% essential_genes_symbol]

