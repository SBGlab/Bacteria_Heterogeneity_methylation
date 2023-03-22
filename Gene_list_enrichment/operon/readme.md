# Methylated gene list operon-based and subsequent functional enrichment analysis

Bedtools intersect of peak methylation bed file list with operons.bed file from [RegulonDB](https://regulondb.ccg.unam.mx/menu/download/full_version/files/11.1/regulonDB11.1_Data_Dist.tar.gz)

Criteria for operon call, 3 or more methylated positions at:
- closest TSS: -1000 to 0 bp (bedtools closest)
- on operon CDS boundries (bedtools intersect)
