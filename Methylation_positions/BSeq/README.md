# Mapping with Bismark

Dataset form SRA [SRP329794](https://www.ncbi.nlm.nih.gov/sra?term=SRP329794), E.coli K12 substr. MG1655 samples used: 
- [SRR15242564](https://www.ncbi.nlm.nih.gov/sra/SRX11548393[accn]): clone 1, WGBS primers (NNNNNN)
- [SRR15242565](https://www.ncbi.nlm.nih.gov/sra/SRX11548394[accn]): clone 1, ABBS primers (NNNNNSuperG)
- [SRR15242568](https://www.ncbi.nlm.nih.gov/sra/SRX11548397[accn]): clone 2, WGBS primers (NNNNNN)
- [SRR15242569](https://www.ncbi.nlm.nih.gov/sra/SRX11548398[accn]): clone 2, ABBS primers (NNNNNSuperG)

Reference genome: [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3)

# Partial 5mC call from clones 1 and 2

A CX_report was extracted with bismak_methylation_extractor. The selection criteria for cytosine positions was a percentaje of methylation between 25 and 75, both included, and a coverage equal or above 25 found in both sample preparations.  
All the hits were found in a **CHG context**.
