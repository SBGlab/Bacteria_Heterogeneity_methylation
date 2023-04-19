# Methylated gene list operon-based and subsequent functional enrichment analysis

Bedtools intersect of peak methylation bed file list with operons.bed file from [RegulonDB](https://regulondb.ccg.unam.mx/menu/download/full_version/files/11.1/regulonDB11.1_Data_Dist.tar.gz)

Criteria for operon call:
- 4 or more methylated positions at CDS (bedtools intersect with operon_cds_sorted.bed)
- 1 or more methylated positions at regulatory area (bedtools intersect with operon_TSS_sorted.bed)

