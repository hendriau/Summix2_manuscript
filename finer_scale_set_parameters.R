setwd("/storage/math/projects/gnomad-public/gnomADv3.1.2/all_chrs/Publication_FST_sims")

source("/home/priceade/2waysims_parallel_sampledat_fn.R")


#read in merged and AF filtered data
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
library(Summix)
library(parallel)
library(LaF)
library(stringr)

data_list <- vector(mode='list', length=100)


for (i in 1:100){
    sampled_dat <- fread(paste0('/storage/math/projects/gnomad-public/gnomADv3.1.2/all_chrs/Bootstrapped_Dataframes/f5&s1', i, '.text.gz'), header=T)
    data_list[[i]] <- sampled_dat
}




###################################### 100 #############################################
sim_path = "results/AFR_genomewide_parallel.f5_s1.scenario"


#Real data sample sizes
N_bantukenya = 100
N_bantusafrica = 100
N_basque = 100
N_beb = 100
N_cambodian = 100
N_cdx = 100
N_ceu = 100
N_chb = 100
N_chs = 100
N_colombian = 100
N_dai = 100
N_daur = 100
N_esn = 100
N_fin = 100
N_french = 100
N_gbr = 100
N_gih = 100
N_gwd = 100
N_han = 100
N_hezhen = 100
N_ibs = 100
N_italian = 100
N_itu = 100
N_japanese = 100
N_jpt = 100
N_karitiana = 100
N_khv = 100
N_lahu = 100
N_lwk = 100
N_mandenka = 100
N_maya = 100
N_miaozu = 100
N_mongola = 100
N_msl = 100
N_naxi = 100
N_orcadian = 100
N_oroqen = 100
N_pima = 100
N_pjl = 100
N_sardinian = 100
N_she = 100
N_stu = 100
N_surui = 100
N_tsi = 100
N_tu = 100
N_tujia = 100
N_tuscan = 100
N_xibo = 100
N_yizu = 100
N_yoruba = 100
N_yri = 100



#define scenario number
scenario = 1
#Create vector of fine-scale ancestries to be simulated
anc_list <- c("bantukenya", "bantusafrica", "basque", "beb", "cambodian", "cdx", "ceu", "chb", "chs", "colombian", "dai", "daur", "esn", "fin", "french", "gbr", "gih", "gwd", "han", "hezhen", "ibs", "italian", "itu", "japanese", "jpt", "karitiana", "khv", "lahu", "lwk", "mandenka", "maya", "miaozu", "mongola", "msl", "naxi", "orcadian", "oroqen", "pima", "pjl", "sardinian", "she", "stu", "tsi", "tu", "tujia", "tuscan", "xibo", "yizu", "yoruba", "yri")

#Create vector of  data sample sizes
sample_N <- c(paste(rep("N_", each = length(anc_list)), anc_list, sep = ""))
sample_counts <- as.data.frame(lapply(sample_N, get))

#Create parameter dataset
col_names <- paste(c("A", "S" , "E"), rep(1:length(anc_list), each = 3), sep = "")
a <- c("bantukenya", "0", "0", "bantusafrica", "0", "0", "basque", "0", "0", "beb", "0", "0", "cambodian", "0", "0", "cdx", "0", "0", "ceu", "0", "0", "chb", ".25", ".25", "chs", "0", "0", "colombian", "0", "0", "dai", "0", "0", "daur", "0", "0", "esn", "0", "0", "fin", "0", "0", "french", "0", "0", "gbr", "0", "0", "gih", "0", "0", "gwd", "0", "0", "han", "0", "0", "hezhen", "0", "0", "ibs", "0", "0", "italian", "0", "0", "itu", "0", "0", "japanese", "0", "0", "jpt", "0", "0", "karitiana", "0", "0", "khv", ".25", ".25", "lahu", "0", "0", "lwk", "0", "0", "mandenka", "0", "0", "maya", "0", "0", "miaozu", "0", "0", "mongola", ".25", ".25", "msl", "0", "0", "naxi", "0", "0", "orcadian", "0", "0", "oroqen", "0", "0", "pima", "0", "0", "pjl", "0", "0", "sardinian", "0", "0", "she", "0", "0", "stu", "0", "0", "tsi", "0", "0", "tu", ".25", ".25", "tujia", "0", "0", "tuscan", "0", "0", "xibo", "0", "0", "yizu", "0", "0", "yoruba", "0", "0", "yri", "0", "0")
params <- data.frame(matrix(a, nrow = 1, ncol = length(col_names), dimnames = list(c(), col_names)))
params[,grep('S|E', names(params))] <- lapply(params[,grep('S|E', names(params))], as.character)
params[,grep('S|E', names(params))] <- lapply(params[,grep('S|E', names(params))], as.numeric)


sim_Fn <- function(x){
    two_way_simulations(anc_list = anc_list, 
			params = params, 
			real_data = x, 
			sample_counts = sample_counts)}


#run simulations in parallel
affinity <- c(1:50, 1:50)
sim_window <- mclapply(data_list, sim_Fn, mc.preschedule = F, affinity.list = affinity, mc.cores = 50L)


saveRDS(sim_window, paste0(sim_path, scenario, ".rds"))









