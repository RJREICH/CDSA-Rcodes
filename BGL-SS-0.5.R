#-------------------------------------------------------------------
# Function for bayesian Group Lasso with Spike and Slab application 
# Prior Probability Inclusion is set to 0.5
#-------------------------------------------------------------------

# Data_simulated_nstudies: Simulated data obtained from Simulated_data R script
# nb_studies: nb of silulated studies
# nb_iter: number of MCMC iterations
# nb_burn: number of iterations to discard at the beginning
# position.categorical.vector: corresponding to the vector "pos" defined in the article 
# levels.categorical.vector: corresponding to the vector "m" defined in the article 


BGL_p_0.5 <- function(Data_simulated_nstudies,nb_studies,nb_iter,nb_burn,position.categorical.vector,levels.categorical.vector){
  
  #load packages needed
  library(doParallel)
  library(foreach)
  
  cores <- detectCores()
  # create clusters
  cl <- makeCluster(cores-3)
  registerDoParallel(cl)
  
  #Import JAGS BSGL-SS
  
  model.jags <- paste(getwd(),"/BGL-SS-0.5.txt", sep="")
  
  
  result = foreach (i=1:nb_studies,.packages = c("R2jags"),.combine = rbind)%dopar% {
    
    set.seed(123)
    
    # The ith simulated data set 
    Data_sim <- Data_simulated_nstudies[[i]]
    
    # JAGS Data 
    
    data = list(m=levels.categorical.vector ,pos = position.categorical.vector,N = nrow( Data_sim),y = Data_sim$y,X = Data_sim[,-1],nb_levels= ncol(Data_sim)-2)
    
    # Parameters to be traced 
    params <- c("beta","tau","prob_inclusion","lambda")
    
    # Run JAGS model
    jags_BGL_SS<- jags(data = data,parameters.to.save = params,n.iter = nb_iter  
                       ,model.file =  model.jags,n.chains = 3,n.burnin = nb_burn)
    
    # The loop output
    cbind(jags_BGL_SS$BUGSoutput$sims.list$beta,jags_BGL_SS$BUGSoutput$sims.list$prob_inclusion, c(jags_BGL_SS$BUGSoutput$summary[,8], rep(NA, length(jags_BGL_SS$BUGSoutput$sims.list$lambda)-length(jags_BGL_SS$BUGSoutput$summary[,8]))))
  }
  #stop clusters
  stopCluster(cl)
  list(Posterior_Effect_Regression =  result[,1:(ncol(Data_simulated_nstudies[[1]])-1)], 
       Posterior_Inlcusion_Probability = result[,ncol(Data_simulated_nstudies[[1]]):(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-3)],
       Rhat =  result[,(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2):ncol(result)])
}


# Application of the  BGL_p_0.5 function with the following features 

# Function inputs
levels.categorical.vector = c(rep(7,7), rep(6,6), rep(3,3),rep(9,9),rep(7,7),rep(4,4))
position.categorical.vector =  c(rep(1,7),rep(8,6), rep(14,3),rep(17,9),rep(26,7), rep(33,4))
load(file = "Data_sim_100obs_100studies.Rdata")

T1.BGL_p_0.5 <- Sys.time()
Result_BGL_p_0.5 = BGL_p_0.5(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 30000,nb_burn = 10000,position.categorical.vector = position.categorical.vector,levels.categorical.vector = levels.categorical.vector)
T2.BGL_p_0.5 <- Sys.time()

# save results 

Time.BGL_p_0.5 <- T2.BGL_p_0.5-T1.BGL_p_0.5
save(Time.BGL_p_0.5, file = "Time.BGL_p_0.5.Rdata")
save(Result_BGL_p_0.5, file = "Result_BGL_p_0.5.Rdata")


#Results of Table 2 in the article: Mean Posterior Inclusion Probability based on the 100 simulated datasets
Mean_Posterior_Inclusion_Probability_BGL_p_0.5 = round(apply(Result_BGL_p_0.5$Posterior_Inlcusion_Probability,2,mean), digits = 2)


