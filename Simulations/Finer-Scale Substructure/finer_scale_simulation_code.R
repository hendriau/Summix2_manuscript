###############################################################################
# Code to create simulated observed and reference group AFs, and evaluate performance of Summix
# Used for finer-scale reference group simulations evaluating summix() and adjAF()
###############################################################################

# Function to generate simulated dataset of observed and reference group AFs, and evaluate 
# how well Summix is able to estimate refernce group mixing proportions
#* @param anc_list character vector of reference groups to be simulated
#* @param params data frame containing reference group to be simulated and proportions to be simulated for each reference group
#* @param real_data data frame containing allele frequencies for each of the reference groups in 'anc_list'
#* @param sample_counts numeric vector of sample sizes to be used when simulating reference groups

#* @return A data frame containing simulated parameters, Summix estimates, and the difference between Summix estimates and parameters


sim_obs_ref_evalSummix <- function(anc_list, params, real_data, sample_counts){
gc()

# Number which calls specific parameters for simulation, increasing
ivalnum = 1

# Number of simulations to run
testnum = 1
# Number of people to simulate
ntot = 10000


# Load in reference data and paramaters
parameters = params
nodecontrol = ceiling(dim(parameters)[1]/20)
tivec = c(1, round(seq(1, dim(parameters)[1], by = dim(parameters)[1]/nodecontrol))[-1], dim(parameters)[1])

# Pulls reference group names from parameter file
AncFrame = data.frame(parameters[,grep('A', names(parameters))])
# Pulls proportion windows from parameter file
SEframe = data.frame(parameters[,grep('S|E', names(parameters))])
# Calculates observed intervals
testint = c(ifelse(ivalnum == 1, tivec[ivalnum], tivec[ivalnum] + 1), tivec[ivalnum+1])

#Create final simulation data frame output
finalframe = data.frame(matrix(vector(), 0, (6 + 3*length(anc_list) + dim(parameters)[2] + 1),
                               dimnames=list(c(), c('P_Num', 'T_Num', 'Seed', names(parameters), 
                                               paste(anc_list, rep("_parameter", each = length(anc_list)), sep = ""),
                                               'Sum_obj', 'Sum_iterations', 'Sum_time', 'Sum_filt',
                                               paste(anc_list, rep("_estimate", each = length(anc_list)), sep = ""),  
                                               paste(anc_list, rep("_accuracy", each = length(anc_list)), sep = "")))), stringsAsFactors=F)



# Simulations
for (m in testint[1]:testint[2]){
  
  # Pull data from reference
  propvec = SEframe[m,]
  ancvec = AncFrame[m,]
  ancnum = dim(AncFrame)[2]
  outframe = finalframe[FALSE,]
  
  for (j in 1:testnum){
    
    # Set Seed
    seed = as.integer(sample.int(1e8, size=1))
    set.seed(seed)

    # Select 100K SNPS from reference data
    refdat = real_data %>% 
        select(CHROM, POS, REF, ALT, paste(rep("AF", each = length(anc_list)), anc_list, sep="")) 
    names(refdat) <- gsub("AF", "", names(refdat))

    
    # Determines proportions of population
    popvecprop = function(propvec){
      prop = numeric(length = (dim(propvec)[2]/2))

      starts = propvec %>% 
        select(starts_with('S'))
      ends = propvec %>% 
        select(starts_with('E'))
      starts <- as.numeric(starts)
      ends <- as.numeric(ends)

      for (i in (1:length(prop))){
        prop[i] = runif(1, min = starts[i], max = ends[i])
      }

      return(prop)
    }


        # Generate vector of proportions of ancestry
        prop = popvecprop(propvec)


        popnum = numeric(ancnum)
        for (i in 1:ancnum){
            popnum[i] = floor((ntot * prop[i]))
          }

        # Initialize matrix for simulated population
        pop_matrix = as.data.frame(matrix(0, nrow = nrow(refdat), ncol = 3))
        # Generate allele counts for population using multinom 
        for (i in 1:ancnum){
          popmatrixadd = t(sapply(refdat[[as.character(ancvec[1,i])]], function(x){x2<-as.numeric(x); rmultinom(1, popnum[i], prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
          pop_matrix = pop_matrix + popmatrixadd
        }
        # MAF threshold to filter by
        MAF_thresh = 0.01

        # Calcluate allele frequencies
        master_frame_gen1 <- data.frame(refdat[,c(1:4)], pop_matrix)
        master_frame_gen1$AF <- (2 * master_frame_gen1[,5] + master_frame_gen1[,6]) / (2 * ntot)
        # Filter by MAF
        master_frame_gen2 <- master_frame_gen1[master_frame_gen1$AF > MAF_thresh & master_frame_gen1$AF < (1-MAF_thresh),]


 ########Simulate new reference data based on N individuals per real populations in HGDP and 1KG###########
        
        #Real data sample sizes
        sample_N <- c(paste(rep("N_", each = length(anc_list)), anc_list, sep = ""))
        sample_counts <- as.data.frame(lapply(sample_N, get))
        
        refsims = as.data.frame(matrix(0, nrow = nrow(refdat), ncol = length(anc_list)))   
        for (i in 1:ancnum){
          refsimcount = t(sapply(refdat[[as.character(ancvec[1,i])]], function(x){x2<-as.numeric(x); rmultinom(1, as.numeric(sample_counts[i]), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
          refsim = (2 * refsimcount[,1] + refsimcount[,2]) / (2 * as.numeric(sample_counts[i]))
          refsims[,i] = refsim
        }

        refdat_sim <- as.data.frame(cbind(refdat[,1], refdat[,2], refdat[,3], refdat[,4], refsims))
        names(refdat_sim) <- c("CHROM", "POS", "REF", "ALT", anc_list)
	refdat_sim_flt <- refdat_sim %>% filter_at(anc_list, any_vars(. > .01 & .< .99))        

        # Pull reference data
        refdatm = refdat_sim_flt %>% 
          select(POS, REF, ALT, all_of(anc_list))
        # Set Observed data
        obsvecm = master_frame_gen2 %>% 
          select(POS, REF, ALT, AF)

        # Merge reference and observed
        mergeframe = merge(refdatm, obsvecm, by = c("POS","REF","ALT")) %>% 
          select(POS, REF, ALT, all_of(anc_list), AF)
        mergeframe[,4:ncol(mergeframe)] <- sapply(mergeframe[,4:ncol(mergeframe)],as.numeric)


    # Run Summix
        R_sum = summix(mergeframe, reference = all_of(anc_list), observed = "AF")

    # Create dataframe with set parameter values for each fine-scale reference ancestry
    ancs <- anc_list
    propinfo <- data.frame(matrix(0, nrow = ncol(as.data.frame(prop)), ncol = length(ancs)))
    for (i in 1:ancnum){
        props = prop[i]
        propinfo[,i] = props
    }
    names(propinfo)<-c(paste((anc_list), rep("_prop", each=length(anc_list)), sep="_"))

    # Create dataframe with Summix estimations for each fine-scale reference ancestry
    R_ests <- R_sum[,5:length(R_sum)]
    Sum_est <- data.frame(matrix(0, nrow = ncol(as.data.frame(prop)), ncol = length(ancs)))
    for (i in 1:ancnum){
        ests = R_ests[i]
        Sum_est[,i] = ests
    }
    names(Sum_est)<-c(paste((anc_list), rep("_est", each=length(anc_list)), sep="_"))

    # Create dataframe with accuracy of Summix estimations (estimations-parameters)
    Sum_acc <- data.frame(matrix(0, nrow = ncol(as.data.frame(prop)), ncol = length(ancs)))
    for (i in 1:ancnum){
        acc = (Sum_est[,i] - propinfo[,i])
        Sum_acc[,i] = acc
    }
    names(Sum_acc)<-c(paste((anc_list), rep("_acc", each=length(anc_list)), sep="_"))

    # Save test info for each simulation
    testinfo = data.frame(P_Num = m, T_Num = j, Seed = seed)
    # Save parameters
    parinfo = parameters[m,]
    # Convert to character
    ancvec[] <- lapply(ancvec, as.character)

    # Bind Data
    outline = cbind(testinfo, parinfo, propinfo, R_sum[1], R_sum[2], R_sum[3], R_sum[4], Sum_est, Sum_acc)
    outframe[j,] = outline
  }
  
  # Bind Data
  finalframe = rbind(finalframe, outframe)
  }

return(finalframe)

}





# Function to generate simulated dataset of observed and reference group AFs for evaluation of adjAF()
#* @param anc_list character vector of reference groups to be simulated
#* @param params data frame containing reference group to be simulated and proportions to be simulated for each reference group
#* @param real_data data frame containing allele frequencies for each of the reference groups in 'anc_list'
#* @param sample_counts numeric vector of sample sizes to be used when simulating reference groups

#* @return A data frame containing the inputted real data and simulated reference and observed group AFs

sim_obs_ref <- function(anc_list, params, real_data, sample_counts){
  
  refdat = real_data %>% 
    select(CHROM, POS, REF, ALT, paste(rep("AF", each = length(anc_list)), anc_list, sep="")) 
  names(refdat) <- gsub("AF", "", names(refdat))
  
  
  # Number which calls specific paramaters for simulation, increasing
  ivalnum = 1
  
  # Number of simulations to run
  testnum = 1
  # Number of people to simulate
  ntot = 10000
  
  
  # Load in reference data and paramaters
  parameters = params
  nodecontrol = ceiling(dim(parameters)[1]/20)
  tivec = c(1, round(seq(1, dim(parameters)[1], by = dim(parameters)[1]/nodecontrol))[-1], dim(parameters)[1])
  
  # Pulls reference group names from parameter file
  AncFrame = data.frame(parameters[,grep('A', names(parameters))])
  # Pulls proportion windows from parameter file
  SEframe = data.frame(parameters[,grep('S|E', names(parameters))])
  # Calculates observed intervals
  testint = c(ifelse(ivalnum == 1, tivec[ivalnum], tivec[ivalnum] + 1), tivec[ivalnum+1])
  
  
  # Simulations
  
  
  # Pull data from reference
  propvec = SEframe[1,]
  ancvec = AncFrame[1,]
  ancnum = dim(AncFrame)[2]
  
  
  
  # Determines proportions of population
  popvecprop = function(propvec){
    prop = numeric(length = (dim(propvec)[2]/2))
    
    starts = propvec %>% 
      select(starts_with('S'))
    ends = propvec %>% 
      select(starts_with('E'))
    starts <- as.numeric(starts)
    ends <- as.numeric(ends)
    
    for (i in (1:length(prop))){
      prop[i] = runif(1, min = starts[i], max = ends[i])
    }
    
    return(prop)
  }
  
  
  # Generate vector of proportions of ancestry
  prop = popvecprop(propvec)
  
  
  popnum = numeric(ancnum)
  for (i in 1:ancnum){
    popnum[i] = floor((ntot * prop[i]))
  }
  
  
  # Set Seed
  seed = as.integer(Sys.time())
  set.seed(seed)
  print(paste0('seed 1', seed))
  
  
  names(refdat) <-  c("CHROM", "POS", "REF", "ALT", anc_list)
  # Initialize matrix for simulated population
  pop_matrix = as.data.frame(matrix(0, nrow = nrow(refdat), ncol = 3))
  # Generate allele counts for population using multinom 
  for (i in 1:ancnum){
    popmatrixadd = t(sapply(refdat[[as.character(ancvec[1,i])]], function(x){x2<-as.numeric(x); rmultinom(1, popnum[i], prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
    pop_matrix = pop_matrix + popmatrixadd
  }
  
  # Calcluate allele frequencies for simulated observed group
  master_frame_gen1 <- data.frame(refdat[,c(1:4)], pop_matrix)
  master_frame_gen1$AF <- (2 * master_frame_gen1[,5] + master_frame_gen1[,6]) / (2 * ntot)
  
  
  ########Simulate new reference data based on N individuals per real populations in HGDP and 1KG###########
  
  #Real data sample sizes
  sample_N <- c(paste(rep("N_", each = length(anc_list)), anc_list, sep = ""))
  sample_counts <- as.data.frame(lapply(sample_N, get))
  
  
  # Set Seed
  seed = as.integer(Sys.time())
  set.seed(seed)
  print(paste0('seed 2', seed))
  
  #simulate allele counts and calculate simulated allele frequencies for each of the reference groups
  refsims = as.data.frame(matrix(0, nrow = nrow(refdat), ncol = length(anc_list)))   
  for (i in 1:ancnum){
    refsimcount = t(sapply(refdat[[as.character(ancvec[1,i])]], function(x){x2<-as.numeric(x); rmultinom(1, as.numeric(sample_counts[i]), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
    refsim = (2 * refsimcount[,1] + refsimcount[,2]) / (2 * as.numeric(sample_counts[i]))
    refsims[,i] = refsim
  }
  
  
  
  mergeframe = cbind.data.frame(real_data, refsims, master_frame_gen1$AF)
  names(mergeframe) = c(names(real_data), paste0('sim_', anc_list), "Simulated_AF")  
  
  
  #return original data, simulated reference group AFs, and simulated observed group AF
  return(mergeframe)
}




