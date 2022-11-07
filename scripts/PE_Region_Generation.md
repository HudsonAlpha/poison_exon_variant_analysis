# PE Region Geneeration
## Summary
This document walks through the R and Bash code used to generate human Poison Exon regions of interest. 

## Input Data
1. `Non_redundant_cass_mm10_summary.txt` converted from [http://zhanglab.c2b2.columbia.edu/data/CortexAS_PNAS_2015/Non_redundant_cass_mm10_summary.xlsx] (Yan et al. 2015)
2. `SFARI-Gene_genes_07-20-2022release_09-06-2022export.csv`, available from [https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes]
3. `./local_data/mouse_mm10_all_exons.tsv` - From UCSC table browser, see comments
4. `./local_data/NCBI_RefSeq_RefSeqAll_human.bed` From UCSC table browser, see comments
5. `./local_data/refseq_mane_select.bed` from UCSC table browser, see comments

## Requirements:
R (w/ tidyverse), we used version 3.6.1
bedtools, we used version 2.28.0
liftOver command line tool from UCSC, available from [https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/]

```R
library(tidyverse)

# all genes from UCSC Genes and Gene Predictions NCBI RefSeq RefSeqAll
all_gene <- read.csv("NCBI_RefSeq_RefSeqAll_human.bed", sep = "\t", header = TRUE)
bed_allgene <- unique(all_gene[,c("chrom","txStart","txEnd","strand","name2")])
write.table(bed_allgene, "human_all_select_exons.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
# save this for gene naming later

# read in AS-NMD file from Yan. Obtained from https://zhanglab.c2b2.columbia.edu/index.php/Cortex_AS, "Annotation of non-redundant cassette exons"
df <- read.csv("Non_redundant_cass_mm10_summary.txt", sep = "\t", header = TRUE)
# filter for just elements with conserved splice sites in hg19
df_conserved <- df[which(df$AG.GU.in.hg19 == 1),]
dim(df_conserved)
[1] 47918    27
# filter for elements where NMD is induced upon inclusion
df_conserved$NMD <- ifelse(grepl("NMD_in", df_conserved$NMD.classification), "1", ifelse(grepl("other", df_conserved$NMD.classification), "-1", "0"))
df_conserved_ASNMD <- df_conserved[which(df_conserved$NMD > 0),]
dim(df_conserved_ASNMD)
[1] 2730   28
#filter output and create an ID field
output <- df_conserved_ASNMD[,c("chrom","chromStart","chromEnd","strand", "UI.len", "exon.len", "DI.len", "NMD.classification", "gene.symbol")]
output$ID <- paste0(output$strand,",",output$UI.len,",",output$exon.len,",",output$DI.len,",",output$NMD.classification, ",", output$gene.symbol)

# filter for those regions that have canonical or semi-canonical start/stop
ASNMD <- output[,c("chrom","chromStart","chromEnd","ID")]
names(ASNMD) <- c("chrom", "chromStart_region", "chromEnd_region", "ID")
# read in mouse exons pulled from UCSC table browser: tb_ccdsGene, table browser query on ccdsGene
exons <- read.csv("mouse_mm10_all_exons.tsv", sep = "\t", header = FALSE)
names(exons) <- c("chrom", "chromStart_exon", "chromEnd_exon", "info", "number", "strand")

# Regions that have a start in canonical mm10
join_chromStart <- left_join(ASNMD, exons, by = c("chrom", "chromStart_region"="chromStart_exon"))
join_chromStart_true <- join_chromStart[which(is.na(join_chromStart$chromEnd_exon) == FALSE),]
# Regions that DO NOT have a start in canonical mm10
join_chromStart_false <- join_chromStart[which(is.na(join_chromStart$chromEnd_exon) == TRUE),]

# Regions that have an end in canonical mm10
join_chromEnd <- left_join(ASNMD, exons, by = c("chrom", "chromEnd_region"="chromEnd_exon"))
join_chromEnd_true <- join_chromEnd[which(is.na(join_chromEnd$chromStart_exon) == FALSE),]
# Regions that DO NOT have an end in canonical mm10
join_chromEnd_false <- join_chromEnd[which(is.na(join_chromEnd$chromStart_exon) == TRUE),]

# filter for unique starts
join_chromStart_unique <- unique(join_chromStart_true[,c("chrom", "chromStart_region", "chromEnd_region", "chromEnd_exon", "strand", "ID")])
# filter for unique ends
join_chromEnd_unique <- unique(join_chromEnd_true[,c("chrom", "chromStart_region", "chromEnd_region", "chromStart_exon", "strand", "ID")])
# Join the two together by the ID
join_total <- full_join(join_chromStart_unique, join_chromEnd_unique, by = "ID")

# calc out the poison exon introns
names(join_total) <- c("chrom", "A", "D", "B", "strand", "ID", "chrom.2", "A.2", "D.2", "C", "strand.2")
join_total_trim <- join_total[,c("chrom", "A", "B", "C", "D", "strand", "ID")]
join_total_trim$ID <- as.character(join_total_trim$ID)
join_total_trim$UI.len <- vapply(strsplit(join_total_trim$ID,","), `[`, 2, FUN.VALUE=character(1))
join_total_trim$DI.len <- vapply(strsplit(join_total_trim$ID,","), `[`, 4, FUN.VALUE=character(1))
join_total_trim$PE.len <- vapply(strsplit(join_total_trim$ID,","), `[`, 3, FUN.VALUE=character(1))
join_total_trim$element_strand <- vapply(strsplit(join_total_trim$ID,","), `[`, 1, FUN.VALUE=character(1))

exon_map_complete <- join_total_trim[complete.cases(join_total_trim), ]

dim(exon_map_complete)
[1] 1742   11
exon_map_complete[,c(2:5,8:10)] = apply(exon_map_complete[, c(2:5,8:10)],2, function(x) as.numeric(x)) #(margin 1 if rows)

# Filter for regions where the intronic math coordinates check out
region_calc_true <- exon_map_complete[which((exon_map_complete$C - exon_map_complete$B) == (exon_map_complete$UI.len + exon_map_complete$DI.len + exon_map_complete$PE.len)),]
# regions that do not add up correctly
region_calc_false <- exon_map_complete[which((exon_map_complete$C - exon_map_complete$B) != (exon_map_complete$UI.len + exon_map_complete$DI.len + exon_map_complete$PE.len)),]

# write out poison exon coordinates based on the coordinate math
region_calc_true$strand <- as.character(region_calc_true$strand)
region_calc_true$PE_region <- ifelse(region_calc_true$strand == "+", paste0(region_calc_true$chrom, ":", region_calc_true$B + region_calc_true$UI.len+1, "-", region_calc_true$C - region_calc_true$DI.len),
        paste0(region_calc_true$chrom, ":",region_calc_true$B + region_calc_true$DI.len+1, "-",region_calc_true$C - region_calc_true$UI.len))

# collect all the regions that don't check out either beginning or end
# non_calc_exons <- full_join(join_chromStart_false, join_chromEnd_false, by = "ID")
# pull just the regions that don't line up and write those to file
non_calc_regions <- anti_join(ASNMD, region_calc_true, by = "ID")

# write region_calc_true to file
write.table(region_calc_true, "PE_determined_regions_NMD_in_only.bed", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# write non_calc_regions to file
write.table(non_calc_regions, "regions_non_calc_overlap_NMD_in_only.bed", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
```

```bash
# take the regions do not have canonical start/stop and extend their boundaries

awk 'BEGIN {FS=OFS="\t"} { print $1, $2, $2+100, $4","$1","$2","$3 }' regions_non_calc_overlap_NMD_in_only.bed > regions_non_calc_overlap_1ex_NMD_in_only.bed
awk 'BEGIN {FS=OFS="\t"} { print $1, $3-100, $3, $4","$1","$2","$3 }' regions_non_calc_overlap_NMD_in_only.bed > regions_non_calc_overlap_2ex_NMD_in_only.bed

# Intersect those regions with refseq select and get alternative boundaries with the RefSeq Select mm10 scaffold
bedtools intersect -a regions_non_calc_overlap_1ex_NMD_in_only.bed -b mouse_mm10_all_select_exons.tsv -wa -wb > regions_non_calc_overlap_1ex_NMD_in_only_intersect.bed
bedtools intersect -a regions_non_calc_overlap_2ex_NMD_in_only.bed -b mouse_mm10_all_select_exons.tsv -wa -wb > regions_non_calc_overlap_2ex_NMD_in_only_intersect.bed
```
```R
# read in the alternative intersected regions
alt_reg1 <- read.csv("regions_non_calc_overlap_1ex_NMD_in_only_intersect.bed", sep = "\t", header = FALSE)
alt_reg2 <- read.csv("regions_non_calc_overlap_2ex_NMD_in_only_intersect.bed", sep = "\t", header = FALSE)
#label columns based on regions, and mm10 intersect
names(alt_reg1) <- c("chrom_bad", "chromStart_bad", "chromEnd_bad", "ID", "chromExon1", "chromStartExon1", "B", "info", "number", "strand1")
names(alt_reg2) <- c("chrom_bad", "chromStart_bad", "chromEnd_bad", "ID", "chromExon2", "C", "chromEndExon2", "info", "number", "strand2")

# map the poison exon coordinates for the regions now with the alternative boundaries
alt_reg1$ID <- as.character(alt_reg1$ID)
alt_reg1$UI.len <- vapply(strsplit(alt_reg1$ID,","), `[`, 2, FUN.VALUE=character(1))
alt_reg1$DI.len <- vapply(strsplit(alt_reg1$ID,","), `[`, 4, FUN.VALUE=character(1))
alt_reg1$PE.len <- vapply(strsplit(alt_reg1$ID,","), `[`, 3, FUN.VALUE=character(1))
alt_reg1$element_strand <- vapply(strsplit(alt_reg1$ID,","), `[`, 1, FUN.VALUE=character(1))

alt_reg1$chrom <- vapply(strsplit(alt_reg1$ID,","), `[`, 7, FUN.VALUE=character(1))
alt_reg1$A <- vapply(strsplit(alt_reg1$ID,","), `[`, 8, FUN.VALUE=character(1))
alt_reg1$D <- vapply(strsplit(alt_reg1$ID,","), `[`, 9, FUN.VALUE=character(1))

alt_reg2$ID <- as.character(alt_reg2$ID)
alt_reg2$UI.len <- vapply(strsplit(alt_reg2$ID,","), `[`, 2, FUN.VALUE=character(1))
alt_reg2$DI.len <- vapply(strsplit(alt_reg2$ID,","), `[`, 4, FUN.VALUE=character(1))
alt_reg2$PE.len <- vapply(strsplit(alt_reg2$ID,","), `[`, 3, FUN.VALUE=character(1))
alt_reg2$element_strand <- vapply(strsplit(alt_reg2$ID,","), `[`, 1, FUN.VALUE=character(1))

alt_reg2$chrom <- vapply(strsplit(alt_reg2$ID,","), `[`, 7, FUN.VALUE=character(1))
alt_reg2$A <- vapply(strsplit(alt_reg2$ID,","), `[`, 8, FUN.VALUE=character(1))
alt_reg2$D <- vapply(strsplit(alt_reg2$ID,","), `[`, 9, FUN.VALUE=character(1))

join_alt_reg_total <- full_join(alt_reg1, alt_reg2, by = "ID")
exon_map_complete_alt_reg <- join_alt_reg_total[complete.cases(join_alt_reg_total), ]
exon_map_complete_alt_reg[,c(7,11:13,22)] = apply(exon_map_complete_alt_reg[, c(7,11:13,22)],2, function(x) as.numeric(x)) #(margin 1 if rows)

# Filter for regions where the intronic math coordinates check out, now with the alternative boundaries
alt_reg_calc_true <- exon_map_complete_alt_reg[which((exon_map_complete_alt_reg$C - exon_map_complete_alt_reg$B) ==
      (exon_map_complete_alt_reg$UI.len.x + exon_map_complete_alt_reg$DI.len.x + exon_map_complete_alt_reg$PE.len.x)),]
# regions that do not add up correctly, now with the alternative boundaries
alt_reg_calc_false <- exon_map_complete_alt_reg[which((exon_map_complete_alt_reg$C - exon_map_complete_alt_reg$B) !=
      (exon_map_complete_alt_reg$UI.len.x + exon_map_complete_alt_reg$DI.len.x + exon_map_complete_alt_reg$PE.len.x)),]

# write out poison exon coordinates based on the coordinate math, now with the alternative boundaries
alt_reg_calc_true$PE_region <- ifelse(alt_reg_calc_true$element_strand.x == "+", paste0(alt_reg_calc_true$chromExon2, ":", alt_reg_calc_true$B + alt_reg_calc_true$UI.len.x+1, "-", alt_reg_calc_true$C - alt_reg_calc_true$DI.len.x),
        paste0(alt_reg_calc_true$chromExon2, ":", alt_reg_calc_true$B + alt_reg_calc_true$DI.len.x+1, "-", alt_reg_calc_true$C - alt_reg_calc_true$UI.len.x))

# read back in those that passed the first calc
pass1 <- read.csv("PE_determined_regions_NMD_in_only.bed", sep = "\t", header = TRUE)
# rename alt_reg_calc_true columns
pass2 <- alt_reg_calc_true[, c("chromExon1", "A.x", "B", "C", "D.x", "strand1", "ID", "UI.len.x", "DI.len.x", "PE.len.x", "element_strand.x", "PE_region")]
names(pass2) <- c("chrom","A","B","C","D","strand","ID","UI.len","DI.len","PE.len","element_strand","PE_region")

pass2$ID <- gsub(",chr.*", "", pass2$ID)

# join pass1 and pass2
total_pass12 <- rbind(pass1,pass2)

df <- read.csv("Non_redundant_cass_mm10_summary.txt", sep = "\t", header = TRUE)
# filter for just elements with conserved splice sites in hg19
df_conserved <- df[which(df$AG.GU.in.hg19 == 1),]
# filter for elements where NMD is induced upon inclusion
df_conserved$NMD <- ifelse(grepl("NMD_in", df_conserved$NMD.classification), "1", ifelse(grepl("other", df_conserved$NMD.classification), "-1", "0"))
df_conserved_ASNMD <- df_conserved[which(df_conserved$NMD > 0),]
#filter output and create an ID field
output <- df_conserved_ASNMD[,c("chrom","chromStart","chromEnd","strand", "UI.len", "exon.len", "DI.len", "NMD.classification", "gene.symbol")]
output$ID <- paste0(output$strand,",",output$UI.len,",",output$exon.len,",",output$DI.len,",",output$NMD.classification, ",", output$gene.symbol)
# filter for those regions that have canonical or semi-canonical start/stop
ASNMD <- output[,c("chrom","chromStart","chromEnd","ID")]
names(ASNMD) <- c("chrom", "chromStart_region", "chromEnd_region", "ID")
never_pass <- anti_join(ASNMD, total_pass12, by ="ID")
write.table(never_pass, "neverPass_Oct23.bed", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

total_pass12$ID_2 <- gsub(".*NMD","NMD",total_pass12$ID)
total_pass12$ID_2 <- gsub(".*other","other",total_pass12$ID_2)
total_pass12$ID_2 <- gsub(",","_",total_pass12$ID_2)
total_pass12$ASNMD_region <- paste0(total_pass12$chrom, ":", total_pass12$B, "-", total_pass12$C)
total_pass12$E <- gsub(".*:", "", total_pass12$PE_region)
total_pass12$F <- gsub(".*-", "", total_pass12$E)
total_pass12$E <- gsub("-.*", "", total_pass12$E)

#create beds of both PE and ASNMD regions
write.table(total_pass12[,c("chrom", "E","F","ID_2")], "designated_poison_exons_oct23.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(total_pass12[,c("chrom", "B","C","ID_2")], "designated_ASNMD_regions_oct23.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


```bash
#Lifted over into human hg38 using the UCSC liftover tool
liftOver  -minMatch=0.1 designated_poison_exons_oct23.bed Canon_Files/mm10ToHg38.over.chain.gz designated_poison_exons_converted_oct23.bed unMapped_poison_exons_oct23
liftOver  -minMatch=0.1 designated_ASNMD_regions_oct23.bed Canon_Files/mm10ToHg38.over.chain.gz designated_ASNMD_regions_converted_oct23.bed unMapped_ASNMD_regions_oct23

# 3. Ensure the regions lifted over are in the appropriate genes
bedtools intersect -a designated_poison_exons_converted_oct23.bed -b human_all_select_genes.bed -wa -wb > designated_poison_exons_allgene_intersect_oct23.bed
awk 'BEGIN { FS=OFS="\t" } { print $1, $2, $3, $4, $9 }' designated_poison_exons_allgene_intersect_oct23.bed | uniq > designated_poison_exons_labeled_oct23.bed

# remove regions that overlap with any human mane transcript: these would be canonical human exons NCBI RefSeq Select and MANE (nbciRefSeqSelect) dataset table
bedtools intersect -v -a designated_poison_exons_labeled_oct23.bed -b refseq_mane_select.bed | uniq | sort -k1,1 -k2,2n > PE_gene_matching_final_not_in_MANE_oct23.bed
# Create introinc boundaries from NCBI_RefSeq_RefSeqAll introns
bedtools intersect -a PE_gene_matching_final_not_in_MANE_oct23.bed -b NCBI_RefSeq_RefSeqAll_human_introns.bed -wa -wb | uniq | sort -k1,1 -k2,2n > full_region_gene_matching_final_not_in_MANE_oct23.bed
```

```R
library(tidyverse)

gene_matched_final <- read.csv("full_region_gene_matching_final_not_in_MANE_oct23.bed", sep = "\t", header = FALSE)

# read in latest OMIM phenotypes release
phenotypes <- read.csv("Canon_Files/genemap2_SAF.txt", sep = "\t", header = TRUE) #from /cluster/lab/gcooper/hg38/CNV_graphing/cnv_annotation_resources/omim/genemap2.txt
SFARI_Genes <- read.csv("SFARI-Gene_genes_07-20-2022release_09-06-2022export.csv", header = TRUE)

names(gene_matched_final) <- c("chrom_PE", "start_PE", "stop_PE", "ID", "gene", "chrom_region", "start_region", "stop_region", "ID2", "frame", "strand")
gene_matched_final$ID <- gsub("NMD_in.*", "NMD_in", gene_matched_final$ID)
gene_matched_final$ID <- gsub("other.*", "other", gene_matched_final$ID)

phenotypes <- phenotypes %>% mutate_all(na_if,"")
# label regions with OMIM phenotypes
matched_phenotypes <- left_join(gene_matched_final, phenotypes, by = c("gene"="Approved.Gene.Symbol"))
# label regions with SFARI gene scores
matched_SFARI <- left_join(matched_phenotypes,SFARI_Genes, by = c("gene"="gene.symbol"))
# create table for paper
paper_table <- matched_SFARI[,c("chrom_PE", "start_PE", "stop_PE", "ID", "gene", "chrom_region", "start_region", "stop_region", "MIM.Number", "Gene.Symbols", "Phenotypes", "genetic.category", "gene.score", "syndromic")]
# complete alt names for those that have none
names(paper_table)[10] <- "total_gene_name"
write.table(paper_table, "paper_table_draft2_oct23.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
```



```bash
while read stage bam outputfile; do
# ran on LSF using bsub -n 1 -R 'rusage[mem=50000]' -o logfile.txt
samtools depth -a -b Canon_Files/Mus_SCN12389_exons_wPE_1based.bed ${bam} > ${outputfile}
done<mouse_dev_bam_list.tsv
```
