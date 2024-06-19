 
# Pets Pooled Analysis

## Introduction

This repository contains the code and supplementary material for the companion animal pooled analysis study presented in the manuscript "A comprehensive profile of the companion animal gut microbiome integrating reference-based and reference-free methods" by Branck et al. 

## Contents
Supplemental Tables
* Supplemental Table 1: 	Public and private data sources of cat, dog, and human gut metagenomic samples.								
* Supplemental Table 2: 	By-sample taxonomic profiles generated by MetaPhlAn 4 for all companion animal gut metagenomes								
* Supplemental Table 3: 	By-sample taxonomic profiles generated by MetaPhlAn 4 for HMP1-II baseline (time 0 samples) gut metagenomic samples 																	
* Supplemental Table 4: 	By-sample taxonomic profiles generated by MetaPhlAn 4 for Madagascar cohort gut metagenomes								
* Supplemental Table 5: 	"This table provides the number of samples that each SGB is found in each host. Only SGBs that had relative abundance > 10^-5 in at least 3 samples in at least one host are included (n = 2,274). If an SGB did not meet the abundance/prevalence criteria in a given host, considered ""not present"" assigned NA (Methods)."								
* Supplemental Table 6:	Taxonomically classified and unclassified Helicobacteraceae spp. are distributed differently across hosts. Values in the table denote prevalence in each host. SGBs present in at least three samples for at least one host are shown. If an SGB did not meet the abundance/prevalence criteria in a given host, considered "not present" assigned NA (Methods).								
* Supplemental Table 7: Taxonomically classified and unclassified Campylobacter spp. are distributed differently across hosts. Values in the table denote prevalence in each host. SGBs present in at least three samples for at least one host are shown. If an SGB did not meet the abundance/prevalence criteria in a given host, considered "not present" assigned NA (Methods).								
* Supplemental Table 8: Feature-wise testing for differences in dog gut microbiome, at the SGB-level, based on housing (facility vs. household). Output for multivariate linear model ran with MaAsLin 2 ran with dog samples only. Model parameters: fixed term = housing, random term = study, reference = facility.						
* Supplemental Table 9: Feature-wise testing for differences in gut microbiome, at the SGB-level, across host species (cats, dogs, and humans). Multivariate linear models ran with MaAsLin 2 for all samples. Model parameters: fixed term = (host species, housing), random term = study, reference = (host, facility).
- the model was ran three times, each with one of the three hosts noted as the reference, in order to make all pairwise comparisons across hosts. Results were combined, duplicate comparison removed, and FDR-adjusted p-values recalculated (noted by column header: fdr.pval)								
* Supplemental Table 10: Feature-wise testing for differences in ARGs across host species (cats, dogs, and humans). Multivariate linear models ran with MaAsLin 2 for all samples. Model parameters: fixed term = (host species), random term = study, reference = (host).
- the model was ran three times, each with one of the three hosts noted as the reference, in order to make all pairwise comparisons across hosts. Results were combined, duplicate comparison removed, and FDR-adjusted p-values recalculated (noted by column header: fdr.pval)
* Supplemental Table 11: Feature-wise testing for differences in ARGs (organized into antibiotic superclasses) across host species (cats, dogs, and humans). Multivariate linear models ran with MaAsLin 2 for all samples. Model parameters: fixed term = (host species), random term = study, reference = (host).
- the model was ran three times, each with one of the three hosts noted as the reference, in order to make all pairwise comparisons across hosts. Results were combined, duplicate comparison removed, and FDR-adjusted p-values recalculated (noted by column header: fdr.pval)							

Custom scripts:
* 
