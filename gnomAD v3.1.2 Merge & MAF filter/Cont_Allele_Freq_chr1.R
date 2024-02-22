
setwd("/storage/math/projects/gnomad-public/gnomADv3.1.2/chr1")

library(data.table)
chr1 <- fread("chr1AFfilteredallNAremoved.txt.gz")

#create 5 continental frequencies (denoted as "relative frequencies" for pop i,..,K)
chr1$AF_AFR <- (chr1$ANesn * chr1$AFesn + chr1$ANgwd * chr1$AFgwd + chr1$ANlwk * chr1$AFlwk + chr1$ANmsl * chr1$AFmsl + chr1$ANyri * chr1$AFyri + chr1$ANbantukenya * chr1$AFbantukenya + chr1$ANbantusafrica * chr1$AFbantusafrica + chr1$ANmandenka * chr1$AFmandenka + chr1$ANyoruba * chr1$AFyoruba)/ (chr1$ANesn + chr1$ANgwd + chr1$ANlwk + chr1$ANmsl + chr1$ANyri + chr1$ANbantukenya + chr1$ANbantusafrica + chr1$ANmandenka + chr1$ANyoruba)
chr1$AF_IAM <- (chr1$ANcolombian * chr1$AFcolombian + chr1$ANkaritiana * chr1$AFkaritiana + chr1$ANmaya * chr1$AFmaya + chr1$ANpima * chr1$AFpima + chr1$ANsurui * chr1$AFsurui) /  (chr1$ANcolombian + chr1$ANkaritiana + chr1$ANmaya + chr1$ANpima + chr1$ANsurui)
chr1$AF_EAS <- (chr1$ANcdx * chr1$AFcdx + chr1$ANchb * chr1$AFchb + chr1$ANchs * chr1$AFchs + chr1$ANjpt * chr1$AFjpt + chr1$ANkhv * chr1$AFkhv + chr1$ANcambodian * chr1$AFcambodian + chr1$ANdai * chr1$AFdai + chr1$ANdaur * chr1$AFdaur + chr1$ANhan * chr1$AFhan + chr1$ANhezhen * chr1$AFhezhen + chr1$ANjapanese * chr1$AFjapanese + chr1$ANlahu * chr1$AFlahu + chr1$ANmiaozu * chr1$AFmiaozu + chr1$ANmongola * chr1$AFmongola + chr1$ANnaxi * chr1$AFnaxi + chr1$ANoroqen * chr1$AForoqen + chr1$ANshe * chr1$AFshe + chr1$ANtu * chr1$AFtu + chr1$ANtujia * chr1$AFtujia + chr1$ANxibo * chr1$AFxibo + chr1$ANyizu * chr1$AFyizu)/ (chr1$ANcdx + chr1$ANchb + chr1$ANchs + chr1$ANjpt + chr1$ANkhv + chr1$ANcambodian + chr1$ANdai + chr1$ANdaur + chr1$ANhan + chr1$ANhezhen + chr1$ANjapanese + chr1$ANlahu  + chr1$ANmiaozu + chr1$ANmongola + chr1$ANnaxi + chr1$ANoroqen + chr1$ANshe + chr1$ANtu + chr1$ANtujia + chr1$ANxibo + chr1$ANyizu)
chr1$AF_EUR <- (chr1$ANceu * chr1$AFceu + chr1$ANgbr * chr1$AFgbr + chr1$ANibs * chr1$AFibs + chr1$ANtsi * chr1$AFtsi + chr1$ANbasque * chr1$AFbasque + chr1$ANfrench * chr1$AFfrench + chr1$ANitalian * chr1$AFitalian + chr1$ANorcadian * chr1$AForcadian + chr1$ANsardinian * chr1$AFsardinian + chr1$ANtuscan * chr1$AFtuscan) / (chr1$ANceu + chr1$ANgbr + chr1$ANibs + chr1$ANtsi + chr1$ANbasque + chr1$ANfrench + chr1$ANitalian  + chr1$ANorcadian + chr1$ANsardinian + chr1$ANtuscan)
chr1$AF_SAS <- (chr1$ANkalash * chr1$AFkalash + chr1$ANbeb * chr1$AFbeb + chr1$ANgih * chr1$AFgih + chr1$ANitu * chr1$AFitu + chr1$ANpjl * chr1$AFpjl + chr1$ANstu * chr1$AFstu) / (chr1$ANkalash + chr1$ANbeb + chr1$ANgih + chr1$ANitu + chr1$ANpjl + chr1$ANstu)
#filter continental groups for rare variants
chr1flt <- chr1[chr1$AF_AFR>0.01 & chr1$AF_AFR<0.99 | chr1$AF_IAM>0.01 & chr1$AF_IAM<0.99 | chr1$AF_EAS>0.01 & chr1$AF_EAS<0.99| chr1$AF_EUR>0.01 & chr1$AF_EUR<0.99 | chr1$AF_SAS>0.01 & chr1$AF_SAS<0.99,]

#remove NA's from above data frame (Where 0's existed)
chr1flt <- na.omit(chr1flt)


#Create table
write.table(chr1flt, file="chr1AFfilteredallNAremoved_cont_anc.txt", sep="\t", quote = FALSE, row.names=F,col.names=T)
