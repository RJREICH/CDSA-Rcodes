#-------------------------------------------------------------------
# Function for bayesian Group Lasso with Spike and Slab application 
# Uniform distribution for Prior Probability Inclusion
#-------------------------------------------------------------------

# Data_simulated_nstudies: Simulated data obtained from Simulated_data R script
# nb_studies: nb of silulated studies
# nb_iter: number of MCMC iterations
# nb_burn: number of iterations to discard at the beginning
# position.categorical.vector: corresponding to the vector "pos" defined in the article 
# levels.categorical.vector: corresponding to the vector "m" defined in the article 


BGL_p_Unifrom <- function(Data_simulated_nstudies,nb_studies,nb_iter,nb_burn,position.categorical.vector,levels.categorical.vector){
  
  #load packages needed
  library(doParallel)
  library(foreach)
  
  cores <- detectCores()
  # create clusters
  cl <- makeCluster(cores-3)
  registerDoParallel(cl)
  
  #Import JAGS BSGL-SS
  
  model.jags <- paste(getwd(),"/BGL-SS-uniform-p.txt", sep="")
  
  
  result = foreach (i=1:nb_studies,.packages = c("R2jags"),.combine = rbind)%dopar% {
    
  set.seed(123)
   
  # The ith simulated data set   
  Data_sim <- Data_simulated_nstudies[[i]]
    
  # JAGS Data
  
  data = list(m=levels.categorical.vector ,pos = position.categorical.vector,N = nrow( Data_sim),y = Data_sim$y,X = Data_sim[,-1],nb_levels= ncol(Data_sim)-2)
  
  # Parameters to be traced 
  params <- c("beta","tau","prob_inclusion","lambda","p")

  # Run JAGS model
   jags_BGL_SS<- jags(data = data,parameters.to.save = params,n.iter = nb_iter  
                              ,model.file =  model.jags,n.chains = 3,n.burnin = nb_burn)
   # The loop output
   cbind(jags_BGL_SS$BUGSoutput$sims.list$beta,jags_BGL_SS$BUGSoutput$sims.list$prob_inclusion, c(jags_BGL_SS$BUGSoutput$summary[,8], rep(NA, length(jags_BGL_SS$BUGSoutput$sims.list$p)-length(jags_BGL_SS$BUGSoutput$summary[,8]))),jags_BGL_SS$BUGSoutput$sims.list$p)
  }
  #stop clusters
  stopCluster(cl)
  list(Posterior_Effect_Regression =  result[,1:(ncol(Data_simulated_nstudies[[1]])-1)], 
       Posterior_Inlcusion_Probability = result[,ncol(Data_simulated_nstudies[[1]]):(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-3)],
       Posterior_Probability_Inclusion = result[,ncol(result)], 
       Posterior_summary =  result[,((ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2))])
}


# Application of the  BGL_p_Unifrom function with the following features 

# Functions inputs
levels.categorical.vector = c(rep(7,7), rep(6,6), rep(3,3),rep(9,9),rep(7,7),rep(4,4))
position.categorical.vector =  c(rep(1,7),rep(8,6), rep(14,3),rep(17,9),rep(26,7), rep(33,4))
load(file = "Data_sim_100obs_100studies.Rdata")

T1.BGL_p_Uniform <- Sys.time()
Result_BGL_p_Uniform = BGL_p_Unifrom(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 30000,nb_burn = 10000,position.categorical.vector = position.categorical.vector,levels.categorical.vector = levels.categorical.vector)
T2.BGL_p_Uniform <- Sys.time()

# save results 

Time.BGL_p_Uniform <- T2.BGL_p_Uniform-T1.BGL_p_Uniform
save(Time.BGL_p_Uniform, file = "Time.BGL_p_Uniform.Rdata")
save(Result_BGL_p_Uniform, file = "Result_BGL_p_Uniform.Rdata")

#Reuslt of Table 2 in the article: Mean Posterior Inclusion Probability based on the 100 simulated datasets
Mean_Posterior_Inclusion_Probability_p_Unif = round(apply(Result_BGL_p_Uniform$Posterior_Inlcusion_Probability,2,mean), digits = 2)

#Result of Figure 2 in the article

# load package
library(ggplot2)

# Posterior estimation of regression effects by simulated study

# split the posterior regression effects relative to each simulated data set
split_effect_regression_by_study = split(as.data.frame(Result_BGL_p_Uniform$Posterior_Effect_Regression),factor(rep(1:100, each = 3000)))

# compute the posterior median regression effects for each simulated data set
median_effect_regression_by_study = lapply(split_effect_regression_by_study, function (x) lapply(x, median))

# Preparing the data for ggplot()
Median_all <- data.frame (m = unlist(median_effect_regression_by_study), index1 = rep(c(1:37),100), index2 =  rep(c(rep("intercept",1),rep("X1",7),rep("X2",6),rep("X3",3), rep("X4",9), rep("X5",7), rep("X6",4)), 100))

#Data containing the true regression effects from which the data was generated
pointdata <- data.frame(index1 = Median_all$index1, y = rep(c(1,0,1,1,2,2,4,4,0,0,0,0,0,0,0,-2,2,0,0,0,0,0,0,0,0,1,0,1,1,1,1,-2,-2,0,0,0,0), 100))
                        

# Figure 2 in the article
p <- ggplot(Median_all, aes(x=index1, y=unlist(median_effect_regression_by_study), group = factor(index1), color = factor(index2))) + geom_boxplot()+ xlab("")+ylab("Posterior Median Regression Effects")+
     theme(axis.text.x = element_blank()) + theme(legend.title = element_blank())+theme_bw()+theme(axis.title.x=element_blank(),
                                                                                                   axis.text.x=element_blank(),
                                                                                                   axis.ticks.x=element_blank())
plot_Beta_posterior_median <- p+geom_point(data = pointdata,aes(x = index1, y = y),bg = "black", colour = "black", shape = 23, size=0.8)+theme(legend.title=element_blank())

