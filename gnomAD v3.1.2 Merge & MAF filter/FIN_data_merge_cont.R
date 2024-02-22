library(data.table)
setwd("/storage/math/projects/gnomad-public/gnomADv3.1.2/chr1")

#read in merged data
mergeframe <- fread("chr1AFfilteredallNAremoved_FIN_edited.txt.gz", header = TRUE)

#read in merged data with cont AFs and save just cont AFs
cont <- fread("chr1AFfilteredallNAremoved_cont_anc.txt.gz", header = TRUE)
conts <- cont[,c(1:4,577:581)]
dat_merge = merge(mergeframe, conts, by = c("CHROM","POS","REF","ALT"))

write.table(dat_merge, file="chr1AFfilteredallNAremoved_cont_anc_FIN_edited.txt", sep="\t", quote = FALSE, row.names=F,col.names=T)