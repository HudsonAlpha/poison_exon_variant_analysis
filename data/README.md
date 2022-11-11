# Files
1. Supp Table 1: tsv w annotations
2. Supp Table 2: PE bed
3. Supp Table 3: PE-containing intron bed
# Format
## Supplementary Table 1 - Annotated table of poison exon cassettes and their containing introns
Tab-separated text file (`.tsv`) with the following columns
1.	PE_chromosome - genomic coordinates of the poison exon cassette 
2.	PE_start - genomic coordinates of the poison exon cassette
3.	PE_end - genomic coordinates of the poison exon cassette
4.	PE_identifier - unique identifier of the poison exon cassette. Format is `PE_[sequential_number]_of_[total_number]_in_[intron_id]` where `intron_id` is as column #11 
5.	Cassette_Yan_2015_Class - Classification of cassette from Yan et al. 2015
6.	gene_symbol - HGNC gene symbol
7.	alt_gene_list - list of alternative gene names
8.	element_chromosome - genomic coordinates of the RefSeq intron containing the poison exon cassette
9.	element_start - genomic coordinates of the RefSeq intron containing the poison exon cassette
10.	element_end - genomic coordinates of the RefSeq intron containing the poison exon cassette
11.	element_identifier - unique identifier of the RefSeq intron from the UCSC RefSeq RefSeqAll table 
12.	MIM_number - OMIM disease ID
13.	MIM_phenotypes - OMIM phenotypes list
14.	SFARI_genetic_category - From SFARI genes 
15.	SFARI_gene_score - From SFARI genes
16.	SFARI_syndromic - From SFARI genes
## Supplementary Table 2 - BED file of poison exon cassette coordinates (hg38)
1. chrom - Supp Table 1 column 1
2. chromStart - Supp Table 1 column 2
3. chromEnd - Supp Table 1 column 3
4. name - Supp Table 1 column 4
## Supplementary Table 3 - BED file of intronic regions containing one or more poison exons (hg38)
1. chromosome - Supp Table 1 column 8
2. start - Supp Table 1 column 9
3. end - Supp Table 1 column 10
4. name - Supp Table 1 column 11
