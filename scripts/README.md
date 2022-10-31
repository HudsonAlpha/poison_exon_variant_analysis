# Poison Exon

## Summary
Summary here.

This contains example code and things for the paper.

## Local Requirements
Local requirements are bcftools, htslib, bedtools, and R. The versions we used are:
bcftools 1.9
htslib 1.9
bedtools 2.28.0
R 3.6.1 (packages: tidyverse, edgeR, data.table, viridis)

### Anvil data acquisition of NYCKidSeq data
Anvil requirements are bcftools, htslib, and bedtools. The versions we used are:
bcftools 1.9
samtools-1.9
htslib 1.9
bedtools 2.30.0

Example input:
#TODO: what is this and what does it go into?

```
entity:NYCKidSeq_id	vcf_in
GEN19-112P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-112P-1-D1.octopus.legacy.vcf.gz
GEN19-113P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-113P-1-D1.octopus.legacy.vcf.gz
GEN19-114P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-114P-1-D1.octopus.legacy.vcf.gz
GEN19-115P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-115P-1-D1.octopus.legacy.vcf.gz
GEN19-116P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-116P-1-D1.octopus.legacy.vcf.gz
GEN19-117P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-117P-1-D1.octopus.legacy.vcf.gz
GEN19-118P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-118P-1-D1.octopus.legacy.vcf.gz
GEN19-119P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-119P-1-D1.octopus.legacy.vcf.gz
GEN19-120P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-120P-1-D1.octopus.legacy.vcf.gz
GEN19-121P-1-D1.octopus.legacy	gs://fc-secure-749755c3-4dac-45d9-9a02-237822a5b313/NYCKidSeq/GEN19-121P-1-D1.octopus.legacy.vcf.gz
```


Anvil WDL:
[https://portal.firecloud.org/?return=anvil#methods/safelker/bcftools_filter_bedtools_intersect/42]
`./anvil/bcftools_filter_bedtools_intersect.wdl`

Docker image:
https://hub.docker.com/r/safelker/richard_cser_anvil
specifically:
https://hub.docker.com/layers/safelker/richard_cser_anvil/21_Dec_21_update/images/sha256-47c4827d3f1de55b9f24a097bb2b912a5877aa494cd82d69985b44c21b534efb?context=explore

#TODO: figure out what's necessary
ex to manually build locally `docker build -t safelker_richard_cser_anvil:Oct_28 .` required a few modifications to build properly (add -y to some apt-get commands; change a github link to https:// not git://; terra-jupyter-base:1.0.0 doesn't seem to build anymore due to dependency changes but terra-jupyter-base:1.0.11 does)

### Merging, annotating and filtering variants from local VCFs to create a table of high-interest variants.
Make tab-separated input_list.txt, as
```
input_vcf_1.vcf.gz  output_dir/filtered_variants_1.vcf
input_vcf_2.vcf.gz  output_dir/filtered_variants_2.vcf
```
where input_vcf_1 and input_vcf_2 are individual input VCFS to extract variants from.
filtered_variants_1.vcf and filtered_variants_2.vcf are temporary individual outputs.


Uses input data:
#TODO: where are these from?
DDG2P_genes_8_1_2021.csv (from ???)
HGNC_conversion.tsv (from ???)
genemap2.txt (from omim)
regions_of_interest.bed


Needed modifications to run script:
Local location of annotations in `vep_config.ini`



#TODO: recommend we change the input so it's just a single list of VCFs and output directory name, then we can basename the filtered VCF names

Overview
For each proband:
   Filter for depth in proband greater than 10 reads
   Filter for allele balance of proband between 0.2 - 0.8
   Isolate just proband record
   Intersect with regions of interest
   bgzip and index
Merge all proband VCFs
Filter to variants seen twice or less in cohort (bcftools)
VEP annotate (including gnomad, topmed, CADD, GERP)
Filter in R for disease-associated genes, population frequency, CADD, GERP, and proper protein effects (process_filter_label.R)

resulting in:
PE_variants_labeled.tsv manageable output
PE_variants_AD_labeled_full.tsv full vcf info
