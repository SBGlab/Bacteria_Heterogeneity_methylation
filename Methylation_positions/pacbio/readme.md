# PacBio datasets
 
 Raw dataset from:
 -  PacificBiosciences / DevNet [8 plex Ecoli Multiplexed Microbial Assembly](https://github.com/PacificBiosciences/DevNet/wiki/8-plex-Ecoli-Multiplexed-Microbial-Assembly). Reference genome: [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3). Processed with: [SMRT Link v11.0](https://downloads.pacbcloud.com/public/software/installers/smrtlink_11.0.0.146107.zip).

Processed dataset:
- gff file [MethSMRT](http://sysbio.gzzoc.com/download/Prokaryota/Ecoli.gff.tar.gz) from the raw data [SRS674093](https://www.ncbi.nlm.nih.gov/sra/?term=SRS674093)

# Methylation call

See the commands employed at SMRTLink_commands.txt

The methylation call (m6A, m5C, m4C) were filtered by IPDRatio>= 4 & fraction methylation between 25-75% to generate an out.gff
