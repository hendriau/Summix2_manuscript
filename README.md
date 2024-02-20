
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Summix2_manuscript

<!-- badges: start -->
<!-- badges: end -->

In this repository you will find files containing the code used to test
Summix2 (bioconductor link), as described by this manuscript: (bioRxiv
link)…. as well as code used to merge and MAF filter the gnomAD v3.1.2
data for analysis.

Below is an outline of the folders in this repository and their
contents.

<br> <br> <br>

<font size="5">[**gnomAD v3.1.2 Merge & MAF
filter**](https://github.com/hendriau/Summix2_manuscript/tree/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter)
folder: </font>

[**chr1_gnomADv3.1.2_merge.sh**](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/chr1_gnomADv3.1.2_merge.sh)
contains the code used to complete the merge of gnomAD v3.1.2 whole
genome variants with Human Genetic Diversity Project (HGDP) and 1000
Genomes Project (1KG) whole genome variants released in gnomAD v3.1.2 as
separate vcf files. Note: This code was used to merge across chromosome
1; the same code- with only chromosome number modified- was used to
merge across chromosomes 2-22.

To read more about this data, proceed to the [**gnomAD v3.1 blog
post**](https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/)
and then to the [**gnomAD v3.1.2 blog
post**](https://gnomad.broadinstitute.org/news/2021-10-gnomad-v3-1-2-minor-release/),
which describes updates made to the gnomAD v3.1 data release.

[**Filter_Allele_Freq_chr1.R**](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/Filter_Allele_Freq_chr1.R)
contains the code used to filter for a MAF \>= .01 in at least one
finer-scale reference group within the merged gnomAD v3.1.2 and HGDP &
1KG data. Note: This code was used to filter across chromosome 1 SNPs;
the same code- with only chromosome number modified- was used to filter
across chromosomes 2-22 SNPs.
