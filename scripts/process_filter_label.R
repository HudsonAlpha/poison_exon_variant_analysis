library(tidyverse)
library(data.table)

# read in arguments
args <- commandArgs(trailingOnly = TRUE)

# read in output_vepannotated_uniq_var_gene.tsv
proband_df_file <- args[1]
proband_df <- read.csv(proband_df_file, sep = "\t", header = TRUE)
# read in output_AF2_merged.vcf
vcf_file <- args[2]
vcf <- read.csv(vcf_file, header = TRUE, sep = "\t")
# read in latest DECIPHER genes release
DDG2P <- read.csv("DDG2P_genes_8_1_2021.csv", header = TRUE)
# read in HGNC gene name conversion file
HGNC_con <- read.csv("HGNC_conversion.tsv", sep = "\t", header = TRUE)
# read in latest OMIM phenotypes release
OMIM <- read.csv("genemap2.txt", sep = "\t", header = TRUE) #/path/to/omim/genemap2.txt

# change filtering criteria to numeric values and set NA values to zero
proband_df$gnomad3_AC <- as.numeric(as.character(proband_df$gnomad3_AC))
proband_df$gnomad3_AC[is.na(proband_df$gnomad3_AC)] <- 0
proband_df$GerpRS[is.na(proband_df$GerpRS)] <- 0
proband_df$GerpRS <- as.numeric(as.character(proband_df$GerpRS))
proband_df$CADD_PHRED <- as.numeric(as.character(proband_df$CADD_PHRED))
proband_df$topmed_AC <- as.numeric(as.character(proband_df$topmed_AC))
proband_df$topmed_AC[is.na(proband_df$topmed_AC)] <- 0

# filter variants
proband_variants_oi <- proband_df[which(
  + proband_df$CADD_PHRED >= 10 &   # filters on greater than
  + proband_df$GerpRS >= 0 &
  + proband_df$gnomad3_AC < 2 & # filters on LESS than
  + proband_df$topmed_AC < 2),]

# label with alterantive gene symbols
proband_variants_oi_labeled <- left_join(proband_variants_oi, HGNC_con, by = c("HGNC_ID"="HGNC.ID"))
proband_variants_oi_labeled$HGNC.number <- gsub("HGNC:", "", proband_variants_oi_labeled$HGNC_ID)
DDG2P$hgnc.id <- as.character(as.numeric(DDG2P$hgnc.id))
# label with DECHIPER classification
proband_variants_oi_labeled_DDG2P <- inner_join(proband_variants_oi_labeled, DDG2P, by = c("HGNC.number"="hgnc.id"))
# label with OMIM phenotypes
proband_variants_oi_labeled_DDG2P_OMIM <- left_join(proband_variants_oi_labeled_DDG2P, OMIM, by = c("Approved.symbol"="Approved.Gene.Symbol"))
proband_variants_oi_labeled_DDG2P_OMIM_NA <- proband_variants_oi_labeled_DDG2P_OMIM %>% mutate_all(na_if,"")
# filter for variants in genes with a phenotype in OMIM
proband_variants_oi_labeled_phenotyped <- proband_variants_oi_labeled_DDG2P_OMIM_NA[which(is.na(proband_variants_oi_labeled_DDG2P_OMIM_NA$Phenotypes) == FALSE),]
# filter for mutation consequence: those with loss of function, and non-codign variants
proband_variants_oi_labeled_phenotyped$mutation.consequence <- as.character(proband_variants_oi_labeled_phenotyped$mutation.consequence)
proband_variants_oi_labeled_phenotyped_lof <- proband_variants_oi_labeled_phenotyped %>% filter(str_detect(mutation.consequence,'loss of function'))
proband_variants_oi_labeled_phenotyped_lof_allele <- proband_variants_oi_labeled_phenotyped_lof[proband_variants_oi_labeled_phenotyped_lof$allelic.requirement %in% c("monoallelic","x-linked dominant","hemizygous"),]
proband_variants_oi_labeled_phenotyped_lof_collapse <- data.table(proband_variants_oi_labeled_phenotyped_lof_allele)
proband_variants_oi_labeled_phenotyped_lof_collapse$severity <- ifelse(grepl("NMD_transcript_variant|intron_variant|non_coding_transcript_variant|non_coding_transcript_exon_variant", proband_variants_oi_labeled_phenotyped_lof_collapse$Consequence, ignore.case = T), "1",
  ifelse(grepl("missense_variant|coding_sequence_variant|frameshift_variant|inframe_deletion|inframe_insertion|stop_lost", proband_variants_oi_labeled_phenotyped_lof_collapse$Consequence, ignore.case = T), "2",
     ifelse(grepl("splice_region_variant|splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained", proband_variants_oi_labeled_phenotyped_lof_collapse$Consequence, ignore.case = T), "3", "4")))
proband_variants_oi_labeled_phenotyped_lof_collapse$severity <- as.numeric(as.character(proband_variants_oi_labeled_phenotyped_lof_collapse$severity))
proband_variants_oi_labeled_phenotyped_lof_severity <- proband_variants_oi_labeled_phenotyped_lof_collapse[which(proband_variants_oi_labeled_phenotyped_lof_collapse$severity == 1),]

# prepare output
proband_variants_oi_output_vcf <- proband_variants_oi_labeled_phenotyped_lof_severity %>% separate(X.Uploaded_variation, c("X.CHROM", "POS", "REF", "ALT"))
info_output <- proband_variants_oi_output_vcf[, c("X.CHROM", "POS", "REF", "ALT","Approved.symbol")]
info_output$POS <- as.numeric(info_output$POS)

# create a column that labels proband variant
vcf$proband <- apply(vcf=="0/1", 1, FUN= function(x) toString(names(x)[x]))
vcf$POS <- as.numeric(vcf$POS)
# label the filtered variants for proband
proband_pos_of_int_vcf <- inner_join(info_output, vcf, by = c("X.CHROM", "POS", "REF", "ALT"))
proband_out <- proband_pos_of_int_vcf %>% mutate_all(na_if,"")
proband_out_only <- proband_out[which(is.na(proband_out$proband) == FALSE),]

# create manageable output
proband_out_trim <- proband_out_only[,c("X.CHROM", "POS", "REF", "ALT", "Approved.symbol", "proband")]
write.table(proband_out_trim, "PE_variants_labeled.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# write the whole vcf
write.table(proband_out_only, "PE_variants_AD_labeled_full.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
