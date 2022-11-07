# Poison Exon

## TODO
2. stephanie add supp table 1, PE bed, ASNMD bed, update `./data/README.md` with column descriptions
3. stephanie/james Clean up language in readmes
6. james Quick test: Docker container, filter_and_annotate.sh, process_filter_label.R
7. stephanie/james Update PE_Region_Generation.md with code that adds new columns to supp table 1 and has correct filename.
8. James update paper methods text from this as necessary.

## Summary
Summary here.

This contains example code and things for the paper.

## Data

Supplementary data for Felker et al.
1. Supplementary Table 1 (PE regions of interest and annotations)
2. Poison Exons bed file (correspons to columns X and X of Supplementary Table 1)
3. AS-NMD regions of interest (correspons to columns X and X of Supplementary Table 1)

For column descriptions, see `data/README.md`

## Scripts
Supplementary methods for Felker et al.
Contains local code (R, Bash) and AnVIL code (Docker, WDL) used to generate PE regions of interest, extract VCF variants based on these regions, annotate, and filter to variants of interest for manual review.

For more detail, see `scripts/README.md`

## Citation

Please cite our magnificent publication.

## Contact

First author: Stephanie Felker, sfelker@hudsonalpha.org  
PI: Dr. Greg Cooper, gcooper@hudsonalpha.org  
Code Maintainer: James Lawlor, jlawlor@hudsonalpha.org  
