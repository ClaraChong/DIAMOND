# DIAMOND: 16S amplicon and shotgun metagenomic analysis 

Statistical analyses conducted using R (version 3.5.0 and 3.6.1)

PERMANOVA (adonis function in vegan, version 2.5.6, 10,000 permutations) was used to quantify the contributions of covariates to the observed variance in microbial beta diversities

MaAsLin2 was used to idetify associations between individual microbial taxa, that were present in ≥ 10% of the samples

Wilcoxon unpaired (rstatix, version 0.3.0) and Kruskal–Wallis tests were used to compare two and > 2 independent groups, respectively. 

Chi-square test of independence and Fisher's exact test were used to determine the association between two categorical variables

All q-values (FDR corrected p-values) were corrected for multiple testing using Benjamini-Hochberg procedure 
