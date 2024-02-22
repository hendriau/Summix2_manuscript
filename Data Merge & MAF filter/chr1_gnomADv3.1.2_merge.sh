cd /storage/math/projects/gnomad-public/gnomADv3.1.2/chr1


#create files with intersection of gnomad data file and HGDP_1000G data file for ease of merge (isolating REF/ALT allele combos present in both files but pulling from hgdp_1kg data)
bcftools isec -p /storage/math/projects/gnomad-public/gnomADv3.1.2/chr1/ -n=2 -w1 gnomad.genomes.v3.1.2.sites.chr1.PASSonly.recode.vcf.gz gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.chr1.recode.vcf.gz
mv 0000.vcf gnomad.genomes.v3.1.2.sites.chr1.PASSonly.intersect.vcf
bgzip gnomad.genomes.v3.1.2.sites.chr1.PASSonly.intersect.vcf
tabix -p vcf gnomad.genomes.v3.1.2.sites.chr1.PASSonly.intersect.vcf.gz


#create files with intersection of gnomad data file and HGDP_1000G data file for ease of merge (isolating REF/ALT allele combos present in both files but pulling from gnomad data)
bcftools isec -p /storage/math/projects/gnomad-public/gnomADv3.1.2/chr1/ -n=2 -w1 gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.chr1.recode.vcf.gz gnomad.genomes.v3.1.2.sites.chr1.PASSonly.recode.vcf.gz
mv 0000.vcf gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.intersect.chr1.vcf
bgzip gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.intersect.chr1.vcf
tabix -p vcf gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.intersect.chr1.vcf.gz

#merge hgdp_1kg summary data intersect with gnomad data intersect, "none" denotes only records with identical REF and ALT alleles
bcftools merge --merge none gnomad.genomes.v3.1.2.sites.chr1.PASSonly.intersect.vcf.gz  gnomad.genomes.v3.1.2.hgdp_tgp.summarydata.PASSonly.intersect.chr1.vcf.gz > gnomad.genomes.v3.1.2.merge.hgdp_tgp_summarydata.chr1.vcf.gz -Oz
tabix -p vcf gnomad.genomes.v3.1.2.merge.hgdp_tgp_summarydata.chr1.vcf.gz

#extract AC, AN, AF, & nonhomalt INFO columns
vcftools --gzvcf gnomad.genomes.v3.1.2.merge.hgdp_tgp_summarydata.chr1.vcf.gz --get-INFO AC_acb --get-INFO AN_acb --get-INFO AF_acb --get-INFO nhomalt_acb --get-INFO AC_adygei --get-INFO AN_adygei --get-INFO AF_adygei --get-INFO nhomalt_adygei --get-INFO AC_asw --get-INFO AN_asw --get-INFO AF_asw --get-INFO nhomalt_asw --get-INFO AC_balochi --get-INFO AN_balochi --get-INFO AF_balochi --get-INFO nhomalt_balochi --get-INFO AC_bantukenya --get-INFO AN_bantukenya --get-INFO AF_bantukenya --get-INFO nhomalt_bantukenya --get-INFO AC_bantusafrica --get-INFO AN_bantusafrica --get-INFO AF_bantusafrica --get-INFO nhomalt_bantusafrica --get-INFO AC_basque --get-INFO AN_basque --get-INFO AF_basque --get-INFO nhomalt_basque --get-INFO AC_beb --get-INFO AN_beb --get-INFO AF_beb --get-INFO nhomalt_beb --get-INFO AC_bedouin --get-INFO AN_bedouin --get-INFO AF_bedouin --get-INFO nhomalt_bedouin --get-INFO AC_brahui --get-INFO AN_brahui --get-INFO AF_brahui --get-INFO nhomalt_brahui --get-INFO AC_burusho --get-INFO AN_burusho --get-INFO AF_burusho --get-INFO nhomalt_burusho --get-INFO AC_cambodian --get-INFO AN_cambodian --get-INFO AF_cambodian --get-INFO nhomalt_cambodian --get-INFO AC_cdx --get-INFO AN_cdx --get-INFO AF_cdx --get-INFO nhomalt_cdx --get-INFO AC_ceu --get-INFO AN_ceu --get-INFO AF_ceu --get-INFO nhomalt_ceu --get-INFO AC_chb --get-INFO AN_chb --get-INFO AF_chb --get-INFO nhomalt_chb --get-INFO AC_chs --get-INFO AN_chs --get-INFO AF_chs --get-INFO nhomalt_chs --get-INFO AC_clm --get-INFO AN_clm --get-INFO AF_clm --get-INFO nhomalt_clm --get-INFO AC_colombian --get-INFO AN_colombian --get-INFO AF_colombian --get-INFO nhomalt_colombian --get-INFO AC_dai --get-INFO AN_dai --get-INFO AF_dai --get-INFO nhomalt_dai --get-INFO AC_daur --get-INFO AN_daur --get-INFO AF_daur --get-INFO nhomalt_daur --get-INFO AC_druze --get-INFO AN_druze --get-INFO AF_druze --get-INFO nhomalt_druze --get-INFO AC_esn --get-INFO AN_esn --get-INFO AF_esn --get-INFO nhomalt_esn --get-INFO AC_fin --get-INFO AN_fin --get-INFO AF_fin --get-INFO nhomalt_fin --get-INFO AC_french --get-INFO AN_french --get-INFO AF_french --get-INFO nhomalt_french --get-INFO AC_gbr --get-INFO AN_gbr --get-INFO AF_gbr --get-INFO nhomalt_gbr --get-INFO AC_gih --get-INFO AN_gih --get-INFO AF_gih --get-INFO nhomalt_gih --get-INFO AC_gwd --get-INFO AN_gwd --get-INFO AF_gwd --get-INFO nhomalt_gwd --get-INFO AC_han --get-INFO AN_han --get-INFO AF_han --get-INFO nhomalt_han --get-INFO AC_hazara --get-INFO AN_hazara --get-INFO AF_hazara --get-INFO nhomalt_hazara --get-INFO AC_hezhen --get-INFO AN_hezhen --get-INFO AF_hezhen --get-INFO nhomalt_hezhen --get-INFO AC_ibs --get-INFO AN_ibs --get-INFO AF_ibs --get-INFO nhomalt_ibs --get-INFO AC_italian --get-INFO AN_italian --get-INFO AF_italian --get-INFO nhomalt_italian --get-INFO AC_itu --get-INFO AN_itu --get-INFO AF_itu --get-INFO nhomalt_itu --get-INFO AC_japanese --get-INFO AN_japanese --get-INFO AF_japanese --get-INFO nhomalt_japanese --get-INFO AC_jpt --get-INFO AN_jpt --get-INFO AF_jpt --get-INFO nhomalt_jpt --get-INFO AC_kalash --get-INFO AN_kalash --get-INFO AF_kalash --get-INFO nhomalt_kalash --get-INFO AC_karitiana --get-INFO AN_karitiana --get-INFO AF_karitiana --get-INFO nhomalt_karitiana --get-INFO AC_khv --get-INFO AN_khv --get-INFO AF_khv --get-INFO nhomalt_khv --get-INFO AC_lahu --get-INFO AN_lahu --get-INFO AF_lahu --get-INFO nhomalt_lahu --get-INFO AC_lwk --get-INFO AN_lwk --get-INFO AF_lwk --get-INFO nhomalt_lwk --get-INFO AC_makrani --get-INFO AN_makrani --get-INFO AF_makrani --get-INFO nhomalt_makrani --get-INFO AC_mandenka --get-INFO AN_mandenka --get-INFO AF_mandenka --get-INFO nhomalt_mandenka --get-INFO AC_maya --get-INFO AN_maya --get-INFO AF_maya --get-INFO nhomalt_maya --get-INFO AC_miaozu --get-INFO AN_miaozu --get-INFO AF_miaozu --get-INFO nhomalt_miaozu --get-INFO AC_mongola --get-INFO AN_mongola --get-INFO AF_mongola --get-INFO nhomalt_mongola --get-INFO AC_mozabite --get-INFO AN_mozabite --get-INFO AF_mozabite --get-INFO nhomalt_mozabite --get-INFO AC_msl --get-INFO AN_msl --get-INFO AF_msl --get-INFO nhomalt_msl --get-INFO AC_mxl --get-INFO AN_mxl --get-INFO AF_mxl --get-INFO nhomalt_mxl --get-INFO AC_naxi --get-INFO AN_naxi --get-INFO AF_naxi --get-INFO nhomalt_naxi --get-INFO AC_orcadian --get-INFO AN_orcadian --get-INFO AF_orcadian --get-INFO nhomalt_orcadian --get-INFO AC_oroqen --get-INFO AN_oroqen --get-INFO AF_oroqen --get-INFO nhomalt_oroqen --get-INFO AC_palestinian --get-INFO AN_palestinian --get-INFO AF_palestinian --get-INFO nhomalt_palestinian --get-INFO AC_pathan --get-INFO AN_pathan --get-INFO AF_pathan --get-INFO nhomalt_pathan --get-INFO AC_pel --get-INFO AN_pel --get-INFO AF_pel --get-INFO nhomalt_pel --get-INFO AC_pima --get-INFO AN_pima --get-INFO AF_pima --get-INFO nhomalt_pima --get-INFO AC_pjl --get-INFO AN_pjl --get-INFO AF_pjl --get-INFO nhomalt_pjl --get-INFO AC_pur --get-INFO AN_pur --get-INFO AF_pur --get-INFO nhomalt_pur --get-INFO AC_russian --get-INFO AN_russian --get-INFO AF_russian --get-INFO nhomalt_russian --get-INFO AC_sardinian --get-INFO AN_sardinian --get-INFO AF_sardinian --get-INFO nhomalt_sardinian --get-INFO AC_she --get-INFO AN_she --get-INFO AF_she --get-INFO nhomalt_she --get-INFO AC_sindhi --get-INFO AN_sindhi --get-INFO AF_sindhi --get-INFO nhomalt_sindhi --get-INFO AC_stu --get-INFO AN_stu --get-INFO AF_stu --get-INFO nhomalt_stu --get-INFO AC_surui --get-INFO AN_surui --get-INFO AF_surui --get-INFO nhomalt_surui --get-INFO AC_tsi --get-INFO AN_tsi --get-INFO AF_tsi --get-INFO nhomalt_tsi --get-INFO AC_tu --get-INFO AN_tu --get-INFO AF_tu --get-INFO nhomalt_tu --get-INFO AC_tujia --get-INFO AN_tujia --get-INFO AF_tujia --get-INFO nhomalt_tujia --get-INFO AC_tuscan --get-INFO AN_tuscan --get-INFO AF_tuscan --get-INFO nhomalt_tuscan --get-INFO AC_uygur --get-INFO AN_uygur --get-INFO AF_uygur --get-INFO nhomalt_uygur --get-INFO AC_xibo --get-INFO AN_xibo --get-INFO AF_xibo --get-INFO nhomalt_xibo --get-INFO AC_yakut --get-INFO AN_yakut --get-INFO AF_yakut --get-INFO nhomalt_yakut --get-INFO AC_yizu --get-INFO AN_yizu --get-INFO AF_yizu --get-INFO nhomalt_yizu --get-INFO AC_yoruba --get-INFO AN_yoruba --get-INFO AF_yoruba --get-INFO nhomalt_yoruba --get-INFO AC_yri --get-INFO AN_yri --get-INFO AF_yri --get-INFO nhomalt_yri --get-INFO AC_biakapygmy --get-INFO AN_biakapygmy --get-INFO AF_biakapygmy --get-INFO nhomalt_biakapygmy --get-INFO AC_papuan --get-INFO AN_papuan --get-INFO AF_papuan --get-INFO nhomalt_papuan --get-INFO AC_melanesian --get-INFO AN_melanesian --get-INFO AF_melanesian --get-INFO nhomalt_melanesian --get-INFO AC_mbutipygmy --get-INFO AN_mbutipygmy --get-INFO AF_mbutipygmy --get-INFO nhomalt_mbutipygmy --get-INFO AC_san --get-INFO AN_san --get-INFO AF_san --get-INFO nhomalt_san --get-INFO AC_non_neuro_nfe --get-INFO AN_non_neuro_nfe --get-INFO AF_non_neuro_nfe --get-INFO nhomalt_non_neuro_nfe --get-INFO AC_non_v2 --get-INFO AN_non_v2 --get-INFO AF_non_v2 --get-INFO nhomalt_non_v2 --get-INFO AC_non_v2_mid --get-INFO AN_non_v2_mid --get-INFO AF_non_v2_mid --get-INFO nhomalt_non_v2_mid --get-INFO AC_non_topmed_sas --get-INFO AN_non_topmed_sas --get-INFO AF_non_topmed_sas --get-INFO nhomalt_non_topmed_sas --get-INFO AC_oth --get-INFO AN_oth --get-INFO AF_oth --get-INFO nhomalt_oth --get-INFO AC_non_neuro_fin --get-INFO AN_non_neuro_fin --get-INFO AF_non_neuro_fin --get-INFO nhomalt_non_neuro_fin --get-INFO AC_non_v2_asj --get-INFO AN_non_v2_asj --get-INFO AF_non_v2_asj --get-INFO nhomalt_non_v2_asj --get-INFO AC_controls_and_biobanks_ami --get-INFO AN_controls_and_biobanks_ami --get-INFO AF_controls_and_biobanks_ami --get-INFO nhomalt_controls_and_biobanks_ami --get-INFO AC_non_topmed_eas --get-INFO AN_non_topmed_eas --get-INFO AF_non_topmed_eas --get-INFO nhomalt_non_topmed_eas --get-INFO AC_non_v2_amr --get-INFO AN_non_v2_amr --get-INFO AF_non_v2_amr --get-INFO nhomalt_non_v2_amr --get-INFO AC_non_neuro_sas --get-INFO AN_non_neuro_sas --get-INFO AF_non_neuro_sas --get-INFO nhomalt_non_neuro_sas --get-INFO AC_non_v2_oth --get-INFO AN_non_v2_oth --get-INFO AF_non_v2_oth --get-INFO nhomalt_non_v2_oth --get-INFO AC_ami --get-INFO AN_ami --get-INFO AF_ami --get-INFO nhomalt_ami --get-INFO AC_non_v2_sas --get-INFO AN_non_v2_sas --get-INFO AF_non_v2_sas --get-INFO nhomalt_non_v2_sas --get-INFO AC_sas --get-INFO AN_sas --get-INFO AF_sas --get-INFO nhomalt_sas --get-INFO AC_non_cancer_eas --get-INFO AN_non_cancer_eas --get-INFO AF_non_cancer_eas --get-INFO nhomalt_non_cancer_eas --get-INFO AC_non_v2_ami --get-INFO AN_non_v2_ami --get-INFO AF_non_v2_ami --get-INFO nhomalt_non_v2_ami --get-INFO AC_non_neuro --get-INFO AN_non_neuro --get-INFO AF_non_neuro --get-INFO nhomalt_non_neuro --get-INFO AC_controls_and_biobanks_eas --get-INFO AN_controls_and_biobanks_eas --get-INFO AF_controls_and_biobanks_eas --get-INFO nhomalt_controls_and_biobanks_eas --get-INFO AC_fin --get-INFO AN_fin --get-INFO AF_fin --get-INFO nhomalt_fin --get-INFO AC_controls_and_biobanks_afr --get-INFO AN_controls_and_biobanks_afr --get-INFO AF_controls_and_biobanks_afr --get-INFO nhomalt_controls_and_biobanks_afr --get-INFO AC_non_topmed_mid --get-INFO AN_non_topmed_mid --get-INFO AF_non_topmed_mid --get-INFO nhomalt_non_topmed_mid --get-INFO AC_non_topmed --get-INFO AN_non_topmed --get-INFO AF_non_topmed --get-INFO nhomalt_non_topmed --get-INFO AC_non_topmed_amr --get-INFO AN_non_topmed_amr --get-INFO AF_non_topmed_amr --get-INFO nhomalt_non_topmed_amr --get-INFO AC_controls_and_biobanks_amr --get-INFO AN_controls_and_biobanks_amr --get-INFO AF_controls_and_biobanks_amr --get-INFO nhomalt_controls_and_biobanks_amr --get-INFO AC_non_neuro_mid --get-INFO AN_non_neuro_mid --get-INFO AF_non_neuro_mid --get-INFO nhomalt_non_neuro_mid --get-INFO AC_non_v2_afr --get-INFO AN_non_v2_afr --get-INFO AF_non_v2_afr --get-INFO nhomalt_non_v2_afr --get-INFO AC_non_cancer_afr --get-INFO AN_non_cancer_afr --get-INFO AF_non_cancer_afr --get-INFO nhomalt_non_cancer_afr --get-INFO AC_controls_and_biobanks_fin --get-INFO AN_controls_and_biobanks_fin --get-INFO AF_controls_and_biobanks_fin --get-INFO nhomalt_controls_and_biobanks_fin --get-INFO AC_non_cancer_ami --get-INFO AN_non_cancer_ami --get-INFO AF_non_cancer_ami --get-INFO nhomalt_non_cancer_ami --get-INFO AC_eas --get-INFO AN_eas --get-INFO AF_eas --get-INFO nhomalt_eas --get-INFO AC_non_topmed_nfe --get-INFO AN_non_topmed_nfe --get-INFO AF_non_topmed_nfe --get-INFO nhomalt_non_topmed_nfe --get-INFO AC_amr --get-INFO AN_amr --get-INFO AF_amr --get-INFO nhomalt_amr --get-INFO AC_non_neuro_ami --get-INFO AN_non_neuro_ami --get-INFO AF_non_neuro_ami --get-INFO nhomalt_non_neuro_ami --get-INFO AC_non_cancer_mid --get-INFO AN_non_cancer_mid --get-INFO AF_non_cancer_mid --get-INFO nhomalt_non_cancer_mid --get-INFO AC_afr --get-INFO AN_afr --get-INFO AF_afr --get-INFO nhomalt_afr --get-INFO AC_non_cancer_sas --get-INFO AN_non_cancer_sas --get-INFO AF_non_cancer_sas --get-INFO nhomalt_non_cancer_sas --get-INFO AC_non_topmed_fin --get-INFO AN_non_topmed_fin --get-INFO AF_non_topmed_fin --get-INFO nhomalt_non_topmed_fin --get-INFO AC_controls_and_biobanks_mid --get-INFO AN_controls_and_biobanks_mid --get-INFO AF_controls_and_biobanks_mid --get-INFO nhomalt_controls_and_biobanks_mid --get-INFO AC_controls_and_biobanks_sas --get-INFO AN_controls_and_biobanks_sas --get-INFO AF_controls_and_biobanks_sas --get-INFO nhomalt_controls_and_biobanks_sas --get-INFO AC_non_v2_eas --get-INFO AN_non_v2_eas --get-INFO AF_non_v2_eas --get-INFO nhomalt_non_v2_eas --get-INFO AC_mid --get-INFO AN_mid --get-INFO AF_mid --get-INFO nhomalt_mid --get-INFO AC_non_cancer_nfe --get-INFO AN_non_cancer_nfe --get-INFO AF_non_cancer_nfe --get-INFO nhomalt_non_cancer_nfe --get-INFO AC_non_topmed_asj --get-INFO AN_non_topmed_asj --get-INFO AF_non_topmed_asj --get-INFO nhomalt_non_topmed_asj --get-INFO AC_asj --get-INFO AN_asj --get-INFO AF_asj --get-INFO nhomalt_asj --get-INFO AC_non_topmed_ami --get-INFO AN_non_topmed_ami --get-INFO AF_non_topmed_ami --get-INFO nhomalt_non_topmed_ami --get-INFO AC_non_cancer --get-INFO AN_non_cancer --get-INFO AF_non_cancer --get-INFO nhomalt_non_cancer --get-INFO AC_non_v2_fin --get-INFO AN_non_v2_fin --get-INFO AF_non_v2_fin --get-INFO nhomalt_non_v2_fin --get-INFO AC_non_neuro_oth --get-INFO AN_non_neuro_oth --get-INFO AF_non_neuro_oth --get-INFO nhomalt_non_neuro_oth --get-INFO AC_non_neuro_asj --get-INFO AN_non_neuro_asj --get-INFO AF_non_neuro_asj --get-INFO nhomalt_non_neuro_asj --get-INFO AC_non_topmed_afr --get-INFO AN_non_topmed_afr --get-INFO AF_non_topmed_afr --get-INFO nhomalt_non_topmed_afr --get-INFO AC_non_neuro_eas --get-INFO AN_non_neuro_eas --get-INFO AF_non_neuro_eas --get-INFO nhomalt_non_neuro_eas --get-INFO AC_non_cancer_amr --get-INFO AN_non_cancer_amr --get-INFO AF_non_cancer_amr --get-INFO nhomalt_non_cancer_amr --get-INFO AC_controls_and_biobanks --get-INFO AN_controls_and_biobanks --get-INFO AF_controls_and_biobanks --get-INFO nhomalt_controls_and_biobanks --get-INFO AC_controls_and_biobanks_oth --get-INFO AN_controls_and_biobanks_oth --get-INFO AF_controls_and_biobanks_oth --get-INFO nhomalt_controls_and_biobanks_oth --get-INFO AC_non_cancer_oth --get-INFO AN_non_cancer_oth --get-INFO AF_non_cancer_oth --get-INFO nhomalt_non_cancer_oth --get-INFO AC_non_topmed_oth --get-INFO AN_non_topmed_oth --get-INFO AF_non_topmed_oth --get-INFO nhomalt_non_topmed_oth --get-INFO AC_non_v2_nfe --get-INFO AN_non_v2_nfe --get-INFO AF_non_v2_nfe --get-INFO nhomalt_non_v2_nfe --get-INFO AC_controls_and_biobanks_nfe --get-INFO AN_controls_and_biobanks_nfe --get-INFO AF_controls_and_biobanks_nfe --get-INFO nhomalt_controls_and_biobanks_nfe --get-INFO AC_non_cancer_asj --get-INFO AN_non_cancer_asj --get-INFO AF_non_cancer_asj --get-INFO nhomalt_non_cancer_asj --get-INFO AC_non_neuro_amr --get-INFO AN_non_neuro_amr --get-INFO AF_non_neuro_amr --get-INFO nhomalt_non_neuro_amr --get-INFO AC_non_neuro_afr --get-INFO AN_non_neuro_afr --get-INFO AF_non_neuro_afr --get-INFO nhomalt_non_neuro_afr --get-INFO AC_non_cancer_fin --get-INFO AN_non_cancer_fin --get-INFO AF_non_cancer_fin --get-INFO nhomalt_non_cancer_fin --get-INFO AC_controls_and_biobanks_asj --get-INFO AN_controls_and_biobanks_asj --get-INFO AF_controls_and_biobanks_asj --get-INFO nhomalt_controls_and_biobanks_asj --get-INFO AC_nfe --get-INFO AN_nfe --get-INFO AF_nfe --get-INFO nhomalt_nfe --out infoexpand_chr1 

#zip expanded file
gzip infoexpand_chr1.INFO

mv infoexpand_chr1.INFO.gz gnomad.genomes.v3.1.2.merge.hgdp_tgp_summarydata.infoexpanded.chr1.vcf.gz



########Move on to Separate R Script for AF and missing values filter############

Rscript Filter_Allele_Freq_chr1.R


gzip chr1AFfilteredallNAremoved.txt