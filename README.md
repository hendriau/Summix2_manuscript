
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

[**Data Merge & MAF
filter**](https://github.com/hendriau/Summix2_manuscript/tree/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter)
folder:

[chr1_gnomADv3.1.2_merge.sh](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/chr1_gnomADv3.1.2_merge.sh)
contains the code used to complete the merge of gnomAD v3.1.2 whole
genome variants with Human Genetic Diversity Project (HGDP) and 1000
Genomes Project (1KG) whole genome variants released in gnomAD v3.1.2 as
separate vcf files.

To read more about this data, proceed to the [gnomAD v3.1 blog
post](https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/)
and then to the [gnomAD v3.1.2 blog
post](https://gnomad.broadinstitute.org/news/2021-10-gnomad-v3-1-2-minor-release/),
which describes updates made to the gnomAD v3.1 data release.

[Filter_Allele_Freq_chr1.R](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/Filter_Allele_Freq_chr1.R)
contains the code used to filter for a MAF \>= .01 in at least one
finer-scale reference group within the merged gnomAD v3.1.2 and HGDP &
1KG data.

[Cont_Allele_Freq_chr1.R](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/Cont_Allele_Freq_chr1.R)
contains the code used to create continental-level allele frequencies
via a weighted average of the finer-scale reference groups in each
continental groupings as well as an additional MAF \>= .01 across
continental-level AFs.

<br> gnomAD v3.1.2 Finnish and HGDP Finnish allele frequencies had the
same names and were mis-merged in previous files. The merge of these
columns was redone and files were MAF \>= .01 across finer-scale
reference groups and continental level groups in the following three
files:

[expand_hgdp_1kg_INFO_fin_redo.sh](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/expand_hgdp_1kg_INFO_fin_redo.sh)

[FIN_data_merge.R](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/FIN_data_merge.R)

[FIN_data_merge_cont.R](https://github.com/hendriau/Summix2_manuscript/blob/main/gnomAD%20v3.1.2%20Merge%20%26%20MAF%20filter/FIN_data_merge_cont.R)

<br> <br>

[**Simulations/Finer-Scale
Substructure**](https://github.com/hendriau/Summix2_manuscript/tree/main/Simulations/Finer-Scale%20Substructure)
folder:

[finer_scale_set_parameters.R](https://github.com/hendriau/Summix2_manuscript/blob/main/Simulations/Finer-Scale%20Substructure/finer_scale_set_parameters.R)

[finer_scale_simulation_code](https://github.com/hendriau/Summix2_manuscript/blob/main/Simulations/Finer-Scale%20Substructure/finer_scale_simulation_code.R)

<br> [**Simulations/Local
Substructure**](https://github.com/hendriau/Summix2_manuscript/tree/main/Simulations/Local%20Substructure)
folder:

[local_simulation_code.R](https://github.com/hendriau/Summix2_manuscript/blob/main/Simulations/Local%20Substructure/local_simulation_code.R)

[simulation_code.R](https://github.com/hendriau/Summix2_manuscript/blob/main/Simulations/Local%20Substructure/simulation_code.R)

<br> <br>

[**Data**](https://github.com/hendriau/Summix2_manuscript/tree/main/Data)
folder: Contains gnomAD v3.1.2, 1KG, and HGDP real AF data across 35,000
variants randomly sampled from LD 0.2 filtered chromosome 19. Also
included are AFs for 50 simulated finer-scale references (simulated with
n=200), 5 simulated continental references (simulated with n=500), and
an observed group simulated with 1KG-Pel-like proportions.
