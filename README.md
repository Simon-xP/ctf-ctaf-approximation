# cTF cTAF Approximation

## Description
Research and data analysis done by me, Simon Hampton, as part of a co-op term, which involves approximation of cTAF which is an indicator of limit of detection of multi-cancer early detection machine learning models. The project first takes data from BAM files of both tumor and healthy normal control samples. We combine the reads to make simulated cancer plasma samples, then perform variant calling with a mutation list. From this we extract the number of alternate and reference alleles in each mutation type and calculate cTF using the a Poisson distribution in the Bayesian Inference model. cTAF is calculated by multiplying by median NGS tumor purity and should be an accurate prediction of limit of detection. 

## Difficulties
Faced a few of challenges with combining the reads in particular, as each small mistake meant a day's work being wasted as I had to regenerate the files (it took many hours each time) and get busy with other projects or areas of this project while the files were combining and upsampling/downsampling. 

## Files
The cTF_cTAF_calculator.py file is for the calculations given a VCF file, and the allele_counts.py performs mutation calling. This is for the company's use only and would be difficult to replicate without sample data from sequencing panels.