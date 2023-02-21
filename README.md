# Poison exon annotations improve the yield of clinically relevant variants in genomic diagnostic testing

## Summary
This repository contains the supplementary data and methods for Felker 2023 et al. Of primary interest is the supplementary data in TSV (Supp. Table 1) and BED (Supp Table 2, 3) format. The tables describe hg38 genomic coordinates of intronic regions that may be affected by 
clinically-relevant poison exon inclusion Regions were identified via conservation to mouse cortex RNAseq data (Yan et al. 2015) and annotated with data from OMIM and SFARI genes to determine potential clinical relevance.

## Data

Supplementary data for Felker et al.
1. Supplementary Table 1 (PE-associated intronic regions of interest and annotations)
2. Poison Exons cassettes BED file (hg38)
3. PE-associated intronic regions BED file (hg38) 

For column descriptions, see `data/README.md`

## Scripts
Supplementary methods for Felker et al. 2023, Poison exon annotations improve the yield of clinically relevant variants in genomic diagnostic testing.
Contains local code (R, Bash) and AnVIL code (Docker, WDL) used to generate PE regions of interest, extract VCF variants based on these regions, annotate, and filter to variants of interest for manual review.

For more detail, see `scripts/README.md`

## Citation

https://www.biorxiv.org/content/10.1101/2023.01.12.523654v1

## Contact

First author: Stephanie Felker, sfelker@hudsonalpha.org  
PI: Dr. Greg Cooper, gcooper@hudsonalpha.org  
Code Maintainer: James Lawlor, jlawlor@hudsonalpha.org  

## Changelog

2022-11-11 - initial public release of repository

2023-01-13 - bioRxiv preprint, tagged as Preprint 2023-01-13

2023-02-21 - Update repository with bioRxiv preprint information. Add HGVS notation to Supplementary Table 1 (see `data/Felker2022SuppelementaryTable1_hg38_PE_cassettes_and_elements_annotated.tsv`)
