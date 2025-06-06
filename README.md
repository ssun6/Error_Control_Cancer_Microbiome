
# Correction for spurious taxonomic assignments of k-mer classifiers in low microbial biomass samples using shuffled sequences 

## Transparency and Reproducibility

This repository contains:         

 * Counts files produced by processing the sequences, see folder: `data`
 * Code to run the analysis component of the project, see folder: `R`

## Abstract
**Background**<br />
With the increased use of shotgun metagenome and metatranscriptome sequences in characterizing the microbiome, accurate taxonomic classification of sequencing reads is essential for interpreting microbial community composition and revealing differential microbial signature between groups. K-mer based classifiers such as Kraken2 provide high speed and sensitivity, and are commonly used for low microbial biomass samples. However, their performance can be compromised by specific sources of error without proper parameter settings and incorporation of controls. 
<br />
<br />
**Methods**<br />
In this study, we analyzed six sequencing datasets of human tumor biopsies with Kraken2 and investigated how shared hash values (i.e., identical hash values across different k-mer) and the structure of reference databases can contribute to false positive taxonomic assignments in low biomass samples. 
<br />
<br />
**Results**<br />
We demonstrated that in samples with high non-microbial DNA noise, the classified taxa of Kraken2 in sequencing reads are significantly correlated with that of shuffled sequences using the default setting. These taxa showed a similar distribution as those overrepresented in the hash table construction of the reference database. Incorporation of controls using shuffled reads can separate significant taxa with more robust differences from those more affected by background noise. Although the confidence thresholds needed to minimize noise varied with taxa, a minimum value of 0.2 can also help reduce misclassifications. 
<br />
<br />
**Conclusion**<br />
Our findings highlighted the need for caution when interpreting low-abundance or unexpected taxa in sequencing datasets of low microbial biomass samples. This work contributes to a more comprehensive understanding of the limitations of k-mer based classification tools and provides practical guidance for improving accuracy in microbiome research.
<br />
<br />