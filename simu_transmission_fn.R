###### Functions for simulating transmission ####

# If you find this code useful, please cite:
# Coll et al. Peacock SJ. Definition of a genetic relatedness cutoff to exclude recent transmission of meticillin-resistant Staphylococcus aureus: a genomic epidemiology analysis. Lancet Microbe. 2020 Dec;1(8):e328-e335. doi: 10.1016/S2666-5247(20)30149-X. PMID: 33313577; PMCID: PMC7721685. (https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30149-X/fulltext)

###*** Functions ***############################################################################################################

## Sample, with replacement, from your own distribution
rMydist <- function(n,x) {
  # n = number of samples
  # x = list of probabilities
  sample(seq(1,length(x),1), size = n, prob = x, replace=T)
}

## General model for how distance from transmitted changes over time
## parameters taken from data_fit.R 
gen_mod_day <- function(day, param_general_fit, l_x = 150){
  ### day = number of days from first sample
  ### l_x = maximum number of SNPs
  ### param_general_fit = parameters from fit to data
  
  x <- seq(0,l_x,1)
  aa_g <- param_general_fit[1]
  bb_g <- param_general_fit[2]
  p_aa <- param_general_fit[3]
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}


#### Transmission simulation
simu_runs <- function(timesv, mu, npat, nruns, param_general_fit, t0dist, max_t0 = 0, gen_mod = gen_mod_day){
  ### timesv = time when sample in days from transmission
  ### mu = mutation rate
  ### npat = number of patients
  ### nruns = number of runs
  ### param_general_fit = parameters of generalised function 
  ### t0dist = distribution within source patient at time zero
  ### max_t0 = indicator as to whether to sample from t0dist or to take the maximum number of SNPs. 0 = baseline = sample. If change to 1 then use maximum number
  ### gen_mod = function type for generalised function
  
  ### STORE / PREP
  store <- matrix(0,nruns,npat)
  
  store_limits <- matrix(0,nruns,5)
  
  ### SOURCE PATIENT
  # Distribution of SNPs timesv days after transmission
  dis_source <- gen_mod(timesv,param_general_fit) # Fixed

  for(j in 1:nruns){
    # Sample from SOURCE distribution to give number of SNPs different at this time point in SOURCE patients
    dis_source_from_transm <- rMydist(npat,dis_source)
    
    ### RECIPIENT PATIENT (of transmitted strain, timesv previously)
    # number of mutations likely to have occured since sample
    mut_number <- rpois(npat, mu * timesv)
    r <- 1/mut_number
    
    # RECIPIENT: how many mutations away if assume above number of mutations is 
    # the mean of an exponential distribution (1/mean = rate in this definition of exp)
    # dis_receiver_from_transm <- c()
    # for(i in 1:npat){
    # dis_receiver_from_transm <- c(dis_receiver_from_transm,rexp(1,1/mut_number[i]))}
     
    if(max_t0 == 0){
    dis_receiver_from_transm <- rexp(n=1*length(r), rate=r) + rMydist(npat,t0dist) # sample from time zero distribution 
    }else{dis_receiver_from_transm <- rexp(n=1*length(r), rate=r) + rep(length(t0prob_dist), npat)} # take maximum possible SNP distance
    
    # DISTANCE SOURCE <- TRANS -> RECIPIENT
    mut_dist <- dis_source_from_transm + dis_receiver_from_transm
    
    # STORE
    store[j,] <- c(mut_dist) # store mutation distance
    
    # Calculate cut offs
    h<-hist(mut_dist,breaks = seq(0,300,1),plot=FALSE) 
    h_cum <- cumsum(h$density)
    l1<-length(which(h_cum < 0.7)) # 70% of transmission events captured
    l2<-length(which(h_cum < 0.9)) # 90% of transmission events captured
    l3<-length(which(h_cum < 0.95)) # 95% of transmission events captured
    l4<-length(which(h_cum < 0.99)) # 99% of transmission events captured
    # STORE
    store_limits[j,] <- c(j,l1,l2,l3,l4)
  }
  
  # OUTPUT
  store_limits <- as.data.frame(store_limits)
  colnames(store_limits) <- c("run","70%","90%","95%","99%")
  m_store_limits <- melt(store_limits,id.vars = "run")
  
  # RETURN
  # limits cutoffs / distribution at these limits / example distance
  return(list(store_limits = m_store_limits, store = store))
}
