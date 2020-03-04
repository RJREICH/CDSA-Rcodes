#-------------------------------------------------------------
# function for application of Bayesian Sparse Group Selection
# Prior Probability Selection is set to 0.5
#--------------------------------------------------------------

# Data_simulated_nstudies: Simulated data obtained from Simulated_data R script
# nb_studies: nb of silulated studies
# nb_iter: number of MCMC iterations
# position.categorical.vector: corresponding to the vector "pos" defined in the article 

BSGS_p_0.5 <- function(Data_simulated_nstudies,nb_studies,nb_iter,position.categorical.vector){
  
  #load packages 
  library(doParallel)
  library(foreach)
  
  cores <- detectCores()
  # create clusters
  cl <- makeCluster(cores-3)
  registerDoParallel(cl)
  
  #Import JAGS BSGL-SS
  
   model.jags <- paste(getwd(),"/BSGS-0.5.txt", sep="")
   
   result = foreach (i=1:nb_studies,.packages = c("R2jags"),.combine = rbind)%dopar% {
    
    set.seed(123)
     
    # The ith simulated data set 
    Data_sim <- Data_simulated_nstudies[[i]]
    
    # JAGS Data
    
    data = list(pos = position.categorical.vector,N = nrow(Data_sim),y = Data_sim$y,X = Data_sim[,-1],nb_levels= ncol(Data_sim)-2)
    
    # Parameters to be traced 
    params <- c("beta","tau","ind_var","ind_within_var")
    
    #Initialization of unknown parameters
    
    lm <- lm(data$y~ as.matrix(data$X)-1)
    Beta.init <- as.numeric(coefficients(lm))
    tau.init <- 1/summary(lm)$sigma^2
    
    inits <- list(list(TT = Beta.init[-1] , tau = tau.init),
                  list(TT = runif(1,0.9,1)*Beta.init[-1] , tau = runif(1,0.9,1)*tau.init),
                  list(TT = runif(1,0.9,1)*Beta.init[-1], tau = runif(1,0.9,1)*tau.init))
    
    
    # Run JAGS model
    jags_BSGS<- jags(data = data,parameters.to.save = params,n.iter = nb_iter  
                       ,model.file =  model.jags,n.chains = 3, inits = inits)
    
                   
    # loop output               
    cbind(jags_BSGS$BUGSoutput$sims.list$beta,jags_BSGS$BUGSoutput$sims.list$ind_var,jags_BSGS$BUGSoutput$sims.list$ind_within_var, c(jags_BSGS$BUGSoutput$summary[,8], rep(NA, length(jags_BSGS$BUGSoutput$sims.list$tau)-length(jags_BSGS$BUGSoutput$summary[,8]))), c(jags_BSGS$BUGSoutput$DIC, rep(NA,nrow(jags_BSGS$BUGSoutput$sims.list$beta)-1)))
   }
   # stop clusters
  stopCluster(cl)
  list(Posterior_Effect_Regression =  result[,1:(ncol(Data_simulated_nstudies[[1]])-1)], 
       Posterior_Inlcusion_Probability_cov = result[,ncol(Data_simulated_nstudies[[1]]):(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-3)],
       Posterior_Inlcusion_Probability_within = result[,(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2):(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2+35)], 
       Rhat =  result[,(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2+36)],
       DIC = as.numeric(result[, ncol(result)]))
  }


#-------------------------------------------------------------------
# Application of the BSGS_p_0.5 function with the following features 
#-------------------------------------------------------------------

# Functions inputs
position.categorical.vector =  c(rep(1,7),rep(8,6), rep(14,3),rep(17,9),rep(26,7), rep(33,4))
load(file = "Data_sim_100obs_100studies.Rdata")


T1.BSGS_p_0.5 <- Sys.time()
Result_BSGS_p_0.5 = BSGS_p_0.5(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 50000,position.categorical.vector = position.categorical.vector)
T2.BSGS_p_0.5 <- Sys.time()

#--------------
# save results 
#--------------

Time.BSGS_p_0.5 <- T2.BSGS_p_0.5-T1.BSGS_p_0.5
save(Time.BSGS_p_0.5, file = "Time.BSGS_p_0.5.Rdata")
save(Result_BSGS_p_0.5, file = "Result_BSGS_p_0.5.Rdata")

#---------------------
# plots visualization 
#---------------------

# load packages
library(ggplot2)
library(plyr)

#----- First pannel of Figure 3 in the article -----

# split the MCMC covariates posterior inlcuison probabilities relative to each simulated data set
split_posterior_indicatorcov_by_study = split(as.data.frame(Result_BSGS_p_0.5$Posterior_Inlcusion_Probability_cov),factor(rep(1:100, each = 3000)))
# compute the mean of the posterior inclusion probabilities of covariates
mean_posterior_indicatorcov_by_study = lapply(split_posterior_indicatorcov_by_study, function (x) lapply(x, mean))


posterior_PIP_cov <- data.frame (pip = unlist(mean_posterior_indicatorcov_by_study), index1 = rep(c(1:36),100), index2 =  rep(c(rep("X1",7),rep("X2",6),rep("X3",3), rep("X4",9), rep("X5",7), rep("X6",4)), 100))

plot_PIP_cov <- ggplot(posterior_PIP_cov, aes(x=index1, y=pip, group = factor(index1), color = factor(index2))) + geom_boxplot()+ylab("Posterior Inclusion Probability")+
                theme_bw()+geom_hline(yintercept=0.5, linetype="dashed", color = "red")+theme(axis.text.y = element_text(size = 10), axis.title.x=element_blank(), axis.title.y=element_text(size = 12),axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(legend.title = element_blank())+ ggtitle("prior inclusion probability = 0.5")
  
#Median probability model  (PIP >= 0.5)                                                                                               
apply(Result_BSGS_p_0.5$Posterior_Inlcusion_Probability_cov,2,mean)

# High frequent model throughout MCMC iterations
ind_posterior_cov <- as.data.frame(Result_BSGS_p_0.5$Posterior_Inlcusion_Probability_cov)
rownames(ind_posterior_cov) <- NULL

model_frequency = count(ind_posterior_cov, vars = colnames(ind_posterior_cov))
ind_most_freq_model <- which(model_frequency$freq ==  max(model_frequency$freq))
most_freq_model <- model_frequency[ind_most_freq_model ,]

# Deviance Information Criterion

mean(Result_BSGS_p_0.5$DIC, na.rm = T)





