###############################################################################
# Code to create a simulated dataset to test summix-local
# Will create AMR-like observed sample  with AFR, EUR, IAM continental-like
# substructure and reference data for AFR, EUR, IAM continental groups
###############################################################################

#* Function to generate simulated dataset for local substructure testing
#* @param numSimulations integer, number of simulation replicates to generate
#* @param chrData reference dataframe, ideally contains single chromosome and 
#* must contain columns for AF_AFR, AF_EUR, AF_IAM (references)
#* @param N_afr integer, number of observed homogenous AFR individuals
#* @param N_eur integer, number of observed homogenous EUR individuals
#* @param N_iam integer, number of observed homogenous IAM individuals
#* @param pi_iam double, proportion of IAM ancestry across chromosome
#* @param pi_eur double, proportion of EUR ancestry across chromosome
#* @param pi_afr double, proportion of AFR ancestry across chromosome
#* @param pi_iam_block double, proportion of IAM ancestry in the changed block
#* @param regionLength integer, length of region with different proportions in bp
#* @param N_afr_ref integer, number of reference AFR individuals
#* @param N_eur_ref integer, number of reference EUR individuals
#* @param N_iam_ref integer, number of reference IAM individuals
#* @return list of dataframes containing simulated datasets
doSimulationsFull <- function(numSimulations = 100, chrData, 
                              N_afr = 5000, N_eur = 5000, N_iam = 5000,
                              pi_iam = 0.25, pi_iam_block = NA, 
                              pi_eur = 0.65, pi_afr = 0.10,
                              regionLength = 400000,
                              N_afr_ref = 660, N_eur_ref = 660, 
                              N_iam_ref = 120) {
  simulations <- vector(mode = "list", length = numSimulations)
  
  # for each simulation replicate
  for(rep in 1:numSimulations) {
    # make vectors of the reference AFs
    af_afr1 <- chrData$AF_AFR
    af_eur1 <- chrData$AF_EUR
    af_iam1 <- chrData$AF_IAM
    
    # simulate genotypes and then calculate AFs 
    
    # first do homogenous observed groups
    genos <- t(sapply(af_afr1, function(x) {
      rmultinom(1, N_afr, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_afr <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N_afr))
    
    genos <- t(sapply(af_eur1, function(x) {
      rmultinom(1, N_eur, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_eur <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N_eur))
    
    genos <- t(sapply(af_iam1, function(x) {
      rmultinom(1, N_iam, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_iam <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N_iam))
    
    # second do reference group simulations
    genos <- t(sapply(af_afr1, function(x) {
      rmultinom(1, N_afr_ref, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_afr_ref <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N_afr_ref))
    
    genos <- t(sapply(af_eur1, function(x) {
      rmultinom(1, N_eur_ref, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_eur_ref <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N_eur_ref))
    
    genos <- t(sapply(af_iam1, function(x) {
      rmultinom(1, N_iam_ref, c(x^2, 2*x*(1-x), (1-x)^2))
    }))
    af_sim_iam_ref <- apply(genos, 1, function(x) ((x[1]*2)+x[2])/(2*N_iam_ref))
    
    # final observed AF by a weighted average of homogenous simulations
    af_sim <- af_sim_afr*pi_afr + af_sim_eur*pi_eur + af_sim_iam*pi_iam
    
    sim1 <- chrData
    sim1$SIM <- af_sim
    sim1$AF_AFR_SIM <- af_sim_afr_ref
    sim1$AF_EUR_SIM <- af_sim_eur_ref
    sim1$AF_IAM_SIM <- af_sim_iam_ref
    
    # if simulating a block with different proportions
    if(!is.logical(pi_iam_block)) {
      indices <- which(sim1$POS > 30000000 & 
                         sim1$POS < 30000000+regionLength)
      sim1[indices,]$SIM <- af_sim_afr[indices]*pi_afr + 
        af_sim_eur[indices]*(pi_eur-(pi_iam_block - pi_iam)) +
        af_sim_iam[indices]*(pi_iam_block)
    }
    sim1$CHR <- 22
    simulations[[rep]] <- sim1
  }
  return(simulations)
}