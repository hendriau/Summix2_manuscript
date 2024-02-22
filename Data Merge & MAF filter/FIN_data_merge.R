
######################AF filtered merge data########################3
library(data.table)
setwd("/storage/math/projects/gnomad-public/gnomADv3.1.2/chr1")

mergeframe <- fread("chr1AFfilteredallNAremoved.txt.gz", header = TRUE)

names(mergeframe)[393] <- c("ACfin_gnomad")
names(mergeframe)[394] <- c("ANfin_gnomad")
names(mergeframe)[395] <- c("AFfin_gnomad")
names(mergeframe)[396] <- c("nhomaltfin_gnomad")

fin_dat <- fread("FIN_INFO.chr1.vcf_NEW.gz", header = TRUE)
names(fin_dat)[5] <- c("ANfin_new")

#merge INFO data to merged and filtered data
new_merge = merge(mergeframe, fin_dat, by = c("CHROM","POS","REF","ALT"))

#Assign correct columns to gnomAD data and HGDP&1KG data
new_merge$ANfin <- new_merge$ANfin_new

new_merge$ACfin_gnomad <- new_merge$gnomad_AC_fin
new_merge$ANfin_gnomad <- new_merge$gnomad_AN_fin
new_merge$AFfin_gnomad <- new_merge$gnomad_AF_fin
new_merge$nhomaltfin_gnomad <- new_merge$gnomad_nhomalt_fin
new_merge$AFfin_gnomad <- as.numeric(new_merge$AFfin_gnomad)

#remove additional gnomAD columns
new_merge$ANfin_new <- NULL
new_merge$gnomad_AC_fin <- NULL
new_merge$gnomad_AN_fin <- NULL
new_merge$gnomad_AF_fin <- NULL
new_merge$gnomad_nhomalt_fin <- NULL

write.table(new_merge, file="chr1AFfilteredallNAremoved_FIN_edited.txt", sep="\t", quote = FALSE, row.names=F,col.names=T)
