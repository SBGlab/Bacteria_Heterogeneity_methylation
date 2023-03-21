library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
library(org.EcK12.eg.db)
library(clusterProfiler)
library(stringr)
library(biodbKegg)
operon=read.csv("operon.txt",header = FALSE, sep = "\n")
mylist_split <- lapply(operon[39:nrow(operon),], function(x) strsplit(x, "\t"))
x=t(as.data.frame(mylist_split))
x=x[,c(2,5,6,7)]
x[,4][x[,4] == "forward"] <- '+'
x[,4][x[,4] == "reverse"] <- '-'
x=as.data.frame(x)
x <- data.frame(new_col = 'NC_000913', x)
rownames(x)=NULL
x <- x[, c(1, 3,4,5,2)]
#NC_000913	2169693	2171727	-	gatR
#NC_000913	2068635	2070474	+	yoeA
#NC_000913	390251	394418	+	yaiT
#NC_000913	3651291	3653713	+	yhiS
#NC_000913	380844	382872	-	yaiX
#NC_000913	1465392	1474013	+	ydbA
write.table(x,file="operons.bed",col.names = FALSE,row.names = FALSE,quote = FALSE,sep="\t")


operon_gene=read.csv("OperonSet.txt",sep="\t",comment = "#",header=FALSE)[,c(1,6)]
txdb <- makeTxDbFromGFF("GCF_000005845.2_ASM584v2_genomic.gff")

operon_overlapping_genes=function(operon_intersection,threshold){
  operon=read.csv(file = operon_intersection,header = FALSE, sep = "\t")[,10]
  operon_meth=table(operon)[which(table(operon)>=threshold)]
  operon_names=names(operon_meth)
  genes=unlist(strsplit(operon_gene[(operon_gene[,1] %in% operon_names),2],","))
  txdb <- makeTxDbFromGFF("GCF_000005845.2_ASM584v2_genomic.gff")
  converted_names=na.exclude(AnnotationDbi::select(txdb, keys=genes, keytype="TXNAME",columns="GENEID",multiVals="first")[,2])
  write.table(converted_names,file=paste0(str_remove(operon_intersection,"overlaps_"),"_genelist.txt"),row.names=FALSE, col.names=FALSE, sep="",quote=FALSE)
  return (list(converted_names,genes))
}


genes_0_0_locus=as.vector(operon_overlapping_genes("operon_overlaps_0_0",3)[[1]])
genes_7_7_locus=as.vector(operon_overlapping_genes("operon_overlaps_7-7",3)[[1]])
genes_SRR1536433_locus=as.vector(operon_overlapping_genes("operon_overlaps_SRR1536433",2)[[1]])
genes_Bseq_locus=as.vector(operon_overlapping_genes("operon_overlaps_Bseq",3)[[1]])
genes_0_0_symbol=as.vector(operon_overlapping_genes("operon_overlaps_0_0",3)[[2]])
genes_7_7_symbol=as.vector(operon_overlapping_genes("operon_overlaps_7-7",3)[[2]])
genes_SRR1536433_symbol=as.vector(operon_overlapping_genes("operon_overlaps_SRR1536433",2)[[2]])
genes_Bseq_symbol=as.vector(operon_overlapping_genes("operon_overlaps_Bseq",3)[[2]])

#Gene intersections

all_intersections <- list()
for (i in 2:4) {
  combos <- combn(list(genes_SRR1536433_locus, genes_0_0_locus,genes_Bseq_locus,genes_7_7_locus), i, simplify = FALSE)
  for (j in combos) {
    int <- Reduce(intersect, j)
    #This is the same as doing (for i=4) intersect(unlist(combos[[1]][1]),unlist(combos[[1]][2])) only for two vectors possible
    names_j <- names(list(genes_SRR1536433=genes_SRR1536433_locus, genes_0_0=genes_0_0_locus, genes_Bseq=genes_Bseq_locus, genes_7_7=genes_7_7_locus))[match(j, list(genes_SRR1536433=genes_SRR1536433_locus, genes_0_0=genes_0_0_locus, genes_Bseq=genes_Bseq_locus, genes_7_7=genes_7_7_locus))]
    name <- paste(names_j, collapse = " vs ")
    all_intersections[[name]] <- int
  }
}

all_intersections
write(all_intersections,file="operon_combos.txt")



metabolic_genes=read.csv("conv_38DC2E424FE51674044440701.txt",sep="\t",colClasses = c("NULL","character","NULL","NULL"))
metabolic_genes=as.vector(metabolic_genes)[[1]]
metabolic_genes_locus=na.exclude(AnnotationDbi::select(txdb, keys=metabolic_genes, keytype="TXNAME",columns="GENEID",multiVals="first")[,2])
functional_enrichment=function(gene_symbol,background=metabolic_genes){
  ego <- enrichGO(as.vector(gene_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = background)
  return(ego)
}

functional_0_0=functional_enrichment(genes_0_0_symbol)
barplot(functional_0_0, showCategory=10,title="functional_0_0_operon")
functional_0_0_summary=data.frame(functional_0_0@result)
functional_7_7=functional_enrichment(genes_7_7_symbol)
barplot(functional_7_7, showCategory=10,title="functional_7_7_operon")
functional_7_7_summary=data.frame(functional_7_7@result)
functional_SRR1536433=functional_enrichment(genes_SRR1536433_symbol)
barplot(functional_SRR1536433, showCategory=10,title="functional_SRR1536433_operon")
functional_SRR1536433_summary=data.frame(functional_SRR1536433@result)
functional_Bseq=functional_enrichment(genes_Bseq_symbol)
barplot(functional_Bseq, showCategory=10,title="functional_Bseq_operon")
functional_Bseq_summary=data.frame(functional_Bseq@result)

KEGG_0_0 <- enrichKEGG(as.vector(genes_0_0_locus),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_0_0,showCategory=10,title="KEGG 0_0 operon")
KEGG_7_7 <- enrichKEGG(as.vector(genes_7_7_locus),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_7_7,showCategory=10,title="KEGG 7_7 operon")
KEGG_Bseq <- enrichKEGG(as.vector(genes_Bseq_locus),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_Bseq,showCategory=10,title="KEGG Bseq operon")
KEGG_SRR1536433 <- enrichKEGG(as.vector(genes_SRR1536433_locus),keyType = "kegg",organism = "eco" ,qvalueCutoff = 0.05,use_internal_data = FALSE,universe = metabolic_genes_locus)
barplot(KEGG_SRR1536433,showCategory=10,title="KEGG SRR1536433 operon")


#Genes related with heterogenity in biomass growth and regulated by methylation
expressed_biomass=read.csv("expressed_biomass.txt",header = FALSE, sep = "\n")
expressed_biomass_0_0=expressed_biomass[58:85,]
expressed_biomass_SRR1536433=expressed_biomass[2:31,]
expressed_biomass_7_7=expressed_biomass[33:56,]
expressed_biomass_Bseq=expressed_biomass[87:124,]
expressed_biomass_0_0_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_0_0, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_SRR1536433_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_SRR1536433, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_7_7_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_7_7, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
expressed_biomass_Bseq_symbol=na.exclude(AnnotationDbi::select(txdb, keys=expressed_biomass_Bseq, keytype="GENEID",columns="TXNAME",multiVals="first")[,2])
Bseq_biomass_FE=enrichGO(as.vector(expressed_biomass_Bseq_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(Bseq_biomass_FE, showCategory=10,title="BSeq biomass operon")
expressed_0_0_biomass_FE=enrichGO(as.vector(expressed_biomass_0_0_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(expressed_0_0_biomass_FE, showCategory=10,title="expressed_0_0 biomass operon")
expressed_7_7_biomass_FE=enrichGO(as.vector(expressed_biomass_7_7_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(expressed_7_7_biomass_FE, showCategory=10,title="expressed_7_7 biomass operon")
SRR1536433_biomass_FE=enrichGO(as.vector(expressed_biomass_SRR1536433_symbol),keyType = "SYMBOL",OrgDb = org.EcK12.eg.db,ont = "BP",qvalueCutoff = 0.05,readable = TRUE,universe = metabolic_genes)
barplot(SRR1536433_biomass_FE, showCategory=10,title="SRR1536433 biomass operon")

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
expressed_not_biomass_0_0=expressed_not_biomass[74:94,]
expressed_not_biomass_SRR1536433=expressed_not_biomass[2:34,]
expressed_not_biomass_7_7=expressed_not_biomass[36:72,]
expressed_not_biomass_Bseq=expressed_not_biomass[96:131,]
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
