
# inputs:
#    input_list.txt (two-column text file with input.vcf and output_dir/output.vcf)
#    regions_of_interest.bed
# outputs:
#      output_vepannotated_uniq_var_gene.tsv containing unique variants, gene symbol, and hgnc info
#      output_AF2_merged.vcf.gz a vcf of all samples merged, filtered to low frequency variants (<=2 het variants)




# Filter for depth in proband greater than 10 reads
# Filter for allele balance of proband between 0.2 - 0.8
# Isolate just first VCF sample record (assume: this is proband)
# Intersect with regions of interest
# bgzip and index

while read vcf output_file; do
  #we parallelized the bcftools filter step with LSF, ex: bsub -R rusage[mem=100] -o logfile.txt <bcftools command>
bcftools filter -i 'FORMAT/FT[0]=\"PASS\" && FORMAT/DP[0] > 10 && \
  GT[0]=\"het\" && (FORMAT/AD[0:0])/(FORMAT/DP[0]) < 0.8 && (FORMAT/AD[0:0])/(FORMAT/DP[0]) > 0.2' ${vcf} -Ov |
  bcftools annotate -x INFO,FORMAT | sed '/^##/d'| cut -f 1-10 |
  sed '1s/^/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n/' |
  bedtools intersect -a stdin -b regions_of_interest.bed -header -wa |
  bgzip > ${output_file} && \
  bcftools index ${output_file};
done < input_list.txt

# merge vcfs
bcftools merge output_dir/*vcf.gz -Oz --force-samples -o output_merged.vcf.gz

# filter down those of low frequency in dataset. below filters for heterozygous variants that have a internal frequency of 2 or less
bcftools view -i 'COUNT(GT=\"het\")<3' output_merged.vcf.gz | bgzip > output_AF2_merged.vcf.gz
# save for parental reference later
zcat output_AF2_merged.vcf.gz > output_AF2_merged.vcf
# drop genotypes, split multiallelic sites into biallelic records (-), sort
bcftools view -G output_AF2_merged.vcf.gz | bcftools norm -m - | bcftools sort > output_AF2_merged_nosamples_normalized_sorted.vcf
# trim headers for tabular manipulation
sed '/^##/d' output_AF2_merged_nosamples_normalized_sorted.vcf > output_AF2_merged_nosamples_normalized_sorted_trimmed.vcf
# create ID field to label variant through vep annotation
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3=$1"_"$2"_"$4"_"$5, $4, $5, $6, $7, $8, $9}' output_AF2_merged_nosamples_normalized_sorted_trimmed.vcf > output_AF2_merged_nosamples_normalized_sorted_trimmed_id.vcf
# fix ID column name
sed 's/CHROM_POS_REF_ALT/ID/' output_AF2_merged_nosamples_normalized_sorted_trimmed_id.vcf > output_AF2_merged_nosamples_normalized_sorted_trimmed_idf.vcf
# select just unique
uniq output_AF2_merged_nosamples_normalized_sorted_trimmed_idf.vcf > output_AF2_merged_nosamples_normalized_sorted_trimmed_idf_uniq.vcf
# run VEP
vep --tab --config vep_config.ini -i output_AF2_merged_nosamples_normalized_sorted_trimmed_idf_uniq.vcf -o output_vepannotated.tsv
# remove header for tabular manipulation
sed '/^##/d' output_vepannotated.tsv > output_vepannotated.trim.tsv
# unique variant, gene symbol, hgnc pairing
awk '!seen[$1,$19,$21]++' output_vepannotated.trim.tsv > output_vepannotated_uniq_var_gene.tsv

# Post-annnotation filtering
Rscript process_filter_label.R output_vepannotated_uniq_var_gene.tsv output_AF2_merged.vcf
# outputs: PE_variants_labeled.tsv and PE_variants_AD_labeled_full.tsv
