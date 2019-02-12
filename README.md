# Timing and plotting functions for undifferentiated sarcomas

## Dependencies

The R code presented here contains functions to time whole genome duplications in real time (scale with age of patients) and molecular time, using all mutations and (C>T)pG. It has a few dependencies and required inputs:
+ copy number calls
+ SNV driver calls
+ SNV calls
+indel calls
+ libraries: GenomicRanges, DPClust and dpclust3p (https://github.com/Wedge-Oxford/dpclust_smchet_docker), Rsamtools, Biostrings, BSgenome
+ genome reference
+ clinical annotations
 

## Content

The R code contains functions for plotting and timing, as well as a pipeline which does the following: 
1. load libraries, genome driver genes
2. read in copy number profiles and snv calls
3. infer whole genome duplication status (mode of the major allele)
4. infer multiplicities of mutations (indels and snvs) using dpclust3p functions
5. subset mutations by context (keeps (C>T)pG)
6. plot counts of all mutations vs. (C>T)pG vs. age of the patients for real time calibration
7. time driver mutations relative to WGD
8. plot ploidy vs. fraction LOH
9. time whole genome duplications 
10. plot timing of WGD and driver relative to WGD (using all mutations and only (C>T)pG)
11. scale molecular time to real time using patient age
12. plot real time timing
13. save session and summary data frame
