
cd /storage/math/projects/gnomad-public/gnomADv3.1.2/chr1

#extract AC, AN, AF, & nonhomalt INFO columns
vcftools --gzvcf gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.intersect.chr1.vcf.gz --get-INFO AN_fin --get-INFO gnomad_AC_fin --get-INFO gnomad_AN_fin --get-INFO gnomad_AF_fin --get-INFO gnomad_nhomalt_fin --out infoexpand_chr1 

#zip expanded file
gzip infoexpand_chr1.INFO

mv infoexpand_chr1.INFO.gz FIN_INFO.chr1.vcf_NEW.gz


