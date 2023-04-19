#Methylation motif extraction

1. Motif identification with multiMotifMaker. 
2. Extension of each position in the bed files up to the lenght of each of the identified motifs and extract the corresponding fasta sequence ([mEpigram](https://github.com/Wang-lab-UCSD/mEpigram) bedToFasta.py) so that the motif would fit the sequence exactly. 
3. Repeat step 2 for the filtered bed (partial methylation 25% - 75%) and for the full bed (>25% methylation)
4. [RSAT prokaryotes](http://rsat.sb-roscoff.fr/dna-pattern_form.cgi) parameters: search strand (direct only), substitutions (1), other: default
5. Comparison of motif degeneracy score in full vs filtered (methylation_motifs.R)
