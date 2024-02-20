###############################################################################
# This is code to create a simulated dataset with observed AMR-like continental 
# level substructure and references for all 50 finer-scale and 5 continental
# level references
# Requires reference data to start from in dataframe refdata
###############################################################################

fine_ref_cols <- c("AFbantukenya", "AFbantusafrica", "AFbasque", "AFbeb", "AFcambodian", "AFcdx", "AFceu", "AFchb", "AFchs", "AFcolombian", "AFdai", "AFdaur", "AFesn", "AFfin", "AFfrench", "AFgbr", "AFgih", "AFgwd", "AFhan", "AFhezhen", "AFibs", "AFitalian", "AFitu", "AFjapanese", "AFjpt", "AFkaritiana", "AFkhv", "AFlahu", "AFlwk", "AFmandenka", "AFmaya", "AFmiaozu", "AFmongola", "AFmsl", "AFnaxi", "AForcadian", "AForoqen", "AFpima", "AFpjl", "AFsardinian", "AFshe", "AFstu", "AFtsi", "AFtu", "AFtujia", "AFtuscan", "AFxibo", "AFyizu", "AFyoruba", "AFyri")

group_names <- c("Simulated_AF", reference_cont, reference_fine)
group_refs <- c("AF_AFR", "AF_EUR", "AF_EAS", "AF_SAS", "AF_IAM", fine_ref_cols)
nToSim <- c(10000, n_ref_cont, n_ref_fine)

prop_afr <- .1
prop_eur <- .65
prop_iam <- .25

simDat <- refdata %>% select(CHROM, POS)

for(i in 1:length(group_names)) {
  if (i == 1) {
    genos <- t(sapply(refdata$AF_AFR, function(x) {
      rmultinom(1, nToSim[i], c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_afr <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*nToSim[i]))
    
    genos <- t(sapply(refdata$AF_IAM, function(x) {
      rmultinom(1, nToSim[i], c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_iam <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*nToSim[i]))
    
    genos <- t(sapply(refdata$AF_EUR, function(x) {
      rmultinom(1, nToSim[i], c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_eur <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*nToSim[i]))
    
    af_sim <- af_sim_afr*prop_afr + af_sim_eur*prop_eur + af_sim_iam*prop_iam
  } else {
    #print(paste0("doing ", group_refs[i-1]))
    genos <- t(sapply(refdata[,group_refs[i-1]], function(x) {
      rmultinom(1, nToSim[i], c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*nToSim[i]))
  }
  simDat <- cbind(simDat, af_sim)
  colnames(simDat)[ncol(simDat)] <- group_names[i]
}