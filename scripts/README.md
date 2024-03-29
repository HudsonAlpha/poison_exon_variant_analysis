# Poison Exon Analysis Scripts

## Summary
Scripts and example code for:
1. Poison Exon Region generation
2. Quality-filtering, annotation, and effect-filtering of VCF variants in PE regions
3. Terra / AnVIL quality filtering and selection of VCF variants in PE regions

These scripts are intended as examples and documentation and may require some alteration for local implementation.
They have not been extensively tested in non-local environments.

## Poison Exon Region Generation
See `PE_Region_generation.md`

## Filtering and Annotation of Variants in PE Regions

### Local Requirements
Local requirements are ensembl-vep, bcftools, htslib, bedtools, and R. The versions we used are:
bcftools 1.9  
htslib 1.9  
bedtools 2.28.0  
R 3.6.1 (packages: tidyverse, edgeR, data.table, viridis)  
vep 102 with CADD plugin  

### Local Data Sources & Versions:
#### VEP Data for annotation
1. vep cache: version 102, merged, GRCh38; available from [Ensembl VEP](https://uswest.ensembl.org/info/docs/tools/vep/index.html)
2. hg38 CADDv1.6 hg38 pre-scored whole genome SNVs and gnomad 3.0 Indels, available from <https://cadd.gs.washington.edu/download>
3. TOPMed Freeze 8 variant frequencies, available from <https://bravo.sph.umich.edu/freeze8/hg38/downloads>
4. gnomad genomes v3.1.1 variant frequencies, available from <https://gnomad.broadinstitute.org/downloads>
5. hg38 GERP scores based on 62 mammalian species from UCSC multiz 100-way alignment (generated by Philipp Rentzsch at Berlin Institutes of Health, available from <https://krishna.gs.washington.edu/download/CADD/gerp/gerp_score2_hg38_MAM.bg.gz>)

#### Vep Configuration
See `.vep_config.ini`

#### Additional Data for process_filter_label.R
1. `DDG2P_genes_8_1_2021.csv` - DDG2P Gene List, latest available from:  
<https://www.deciphergenomics.org/redirect?to=https%3A%2F%2Fwww.ebi.ac.uk%2Fgene2phenotype%2Fdownloads%2FDDG2P.csv.gz>
2. `.local_data/HGNC_conversion.tsv` - HGNC id conversion table, from BiomaRt
3. `genemap2.txt` - OMIM gene map, available from <https://www.omim.org/downloads>


### Merging, annotating and filtering variants from local VCFs to create a table of high-interest variants.
Using `filter_and_annotate.sh` and `input_list.txt` we filtered based on variant quality, annotated the variants with VEP (see `vep_config.ini`) and filtered for annotations of interest (see `process_filter_label.R`).

#### Input
`input_list.txt`, a list of VCFs to process (assuming VCFs are single-sample or that the proband is the first sample column)
```
/path/to/sample1.vcf.gz
/path/to/sample2.vcf.gz
/path/to/sample3.vcf.gz
```
#### Overview  
For each proband:  
.....Filter for depth in proband greater than 10 reads  
.....Filter for allele balance of proband between 0.2 - 0.8  
.....Isolate just proband record  
.....Intersect with regions of interest  
.....bgzip and index  
Merge all proband VCFs  
Filter to variants seen twice or less in cohort (bcftools)  
VEP annotate (including gnomad, topmed, CADD, GERP)  
Filter in R for disease-associated genes, population frequency, CADD, GERP, and proper protein effects   (`process_filter_label.R`)


#### Output
`PE_variants_labeled.tsv` brief output table for manual review  
`PE_variants_AD_labeled_full.tsv` full vcf info for manual review

### Acquisition of NYCKidSeq data via AnVIL
We acquired NYCKidSeq whole genome data via AnVIL and developed an AnVIL workflow to extract variants in poison exon regions from each participant VCF.

Merging and downstream analysis were completed locally. Overall:

1. VCFs were copied from the NYCKidseq Anvil bucket into a personal workspace bucket
2. PE Region BED file was uploaded to personal workspace bucket
3. WDL was added as a workflow to the Anvil workflow repository:  
<https://portal.firecloud.org/?return=anvil#methods/safelker/bcftools_filter_bedtools_intersect/42>
4. Docker image was built and registered on DockerHub:  
 <https://hub.docker.com/layers/safelker/richard_cser_anvil/21_Dec_21_update/>
5. Input data list was created based on Anvil template and uploaded to personal workspace bucket.
6. Workflow was instantiated via personal AnVIL workspace
7. Results (quality and PE-region filtered VCFs) were downloaded from personal Anvil workspace
8. Further processing (merging, annotation, and analysis) were completed locally, as above


#### Requirements
Anvil requirements are bcftools, htslib, and bedtools. The versions we used are:
bcftools 1.9  
samtools-1.9  
htslib 1.9  
bedtools 2.30.0  

#### Docker Image
We built a docker image available at  
<https://hub.docker.com/layers/safelker/richard_cser_anvil/21_Dec_21_update/>  

See `./anvil/Dockerfile` for an equivalent buildable dockerfile.
The Dockerfile includes additional tools (such as jupyter) which are
not required for this analysis. These additional items were included to make a general-purpose image for the CSER consoritum AnVIL environment.
Dockerfile developed by Richard Green (UW) and Stephanie Felker (HudonsAlpha) and organized for publishing by James Lawlor (HudsonAlpha).
Note that the Terra base image for analysis was v1.0.0 but updated to v1.0.11 for publication.

#### Workflow
We built a WDL `./anvil/bcftools_filter_bedtools_intersect.wdl` to quality-filter variants and extract variants in poison exon regions.
We registered the workflow with AnVIL here:  
<https://portal.firecloud.org/?return=anvil#methods/safelker/bcftools_filter_bedtools_intersect/42>

#### Input Data Table Example
```
entity:NYCKidSeq_id	vcf_in
nyckidseq_id_1  gs://path-to-personal-bucket/nyckidseq_id_1.vcf.gz
nyckidseq_id_2  gs://path-to-personal-bucket/nyckidseq_id_2.vcf.gz
```
