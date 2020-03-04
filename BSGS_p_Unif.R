#-------------------------------------------------------------
# function for application of Bayesian Sparse Group Selection
# Uniform prior is set to the Prior Probability Selection 
#--------------------------------------------------------------

# Data_simulated_nstudies: Simulated data obtained from Simulated_data R script
# nb_studies: nb of silulated studies
# nb_iter: number of MCMC iterations
# position.categorical.vector: corresponding to the vector "pos" defined in the article 


BSGS_p_unif <- function(Data_simulated_nstudies,nb_studies,nb_iter,position.categorical.vector){
  
  #load packages needed
  library(doParallel)
  library(foreach)
  
  cores <- detectCores()
  # create cluster
  cl <- makeCluster(cores-3)
  registerDoParallel(cl)
  
  #Import JAGS BSGL-SS
  
  model.jags <- paste(getwd(),"/BSGS-unif.txt", sep="")
  
  result = foreach (i=1:nb_studies,.packages = c("R2jags"),.combine = rbind)%dopar% {
    
    set.seed(123)
    
    # The ith simulated data set 
    Data_sim <- Data_simulated_nstudies[[i]]
    
    # JAGS Data
    
    data = list(pos = position.categorical.vector,N = nrow(Data_sim),y = Data_sim$y,X = Data_sim[,-1],nb_levels= ncol(Data_sim)-2)
    
    # Parameters to be traced 
    params <- c("beta","tau","ind_var","ind_within_var","p_var", "p_within")
    
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
   cbind(jags_BSGS$BUGSoutput$sims.list$beta,jags_BSGS$BUGSoutput$sims.list$ind_var,jags_BSGS$BUGSoutput$sims.list$ind_within_var,jags_BSGS$BUGSoutput$sims.list$p_var, jags_BSGS$BUGSoutput$sims.list$p_within,c(jags_BSGS$BUGSoutput$summary[,8], rep(NA, length(jags_BSGS$BUGSoutput$sims.list$tau)-length(jags_BSGS$BUGSoutput$summary[,8]))), c(jags_BSGS$BUGSoutput$DIC, rep(NA,nrow(jags_BSGS$BUGSoutput$sims.list$beta)-1)))
  }
  # stop clusters
  stopCluster(cl)
  list(Posterior_Effect_Regression =  result[,1:(ncol(Data_simulated_nstudies[[1]])-1)], 
       Posterior_Inlcusion_Probability_cov = result[,ncol(Data_simulated_nstudies[[1]]):(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-3)],
       Posterior_Inlcusion_Probability_within = result[,(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2):(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2+35)], 
       Rhat =  result[,(ncol(Data_simulated_nstudies[[1]])+ncol(Data_simulated_nstudies[[1]])-2+36)],
       posterior_probabilty_var_levels = result[,(ncol(result)-2):(ncol(result)-1)],
       DIC = as.numeric(result[, ncol(result)]))
}

#--------------------------------------------------------------------
# Application of the BSGS_p_unif function with the following features 
#--------------------------------------------------------------------

# Functions inputs
position.categorical.vector =  c(rep(1,7),rep(8,6), rep(14,3),rep(17,9),rep(26,7), rep(33,4))
load(file = "Data_sim_100obs_100studies.Rdata")


T1.BSGS_p_unif <- Sys.time()
Result_BSGS_p_unif = BSGS_p_unif(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 50000,position.categorical.vector = position.categorical.vector)
T2.BSGS_p_unif <- Sys.time()

#--------------
# save results 
#--------------

Time.BSGS_p_unif <- T2.BSGS_p_unif-T1.BSGS_p_unif
save(Time.BSGS_p_unif, file = "Time.BSGS_p_unif.Rdata")
save(Result_BSGS_p_unif, file = "Result_BSGS_p_unif.Rdata")

#---------------------
# plots visualization
#---------------------

# load packages 
library(ggplot2)
library(plyr)

# second pannel of Figure 3 in the article

# split the MCMC covariates posterior inlcuison probabilities relative to each simulated data set
split_posterior_indicatorcov_by_study = split(as.data.frame(Result_BSGS_p_unif$Posterior_Inlcusion_Probability_cov),factor(rep(1:100, each = 3000)))
# compute the mean of the posterior inclusion probabilities of covariates
mean_posterior_indicatorcov_by_study = lapply(split_posterior_indicatorcov_by_study, function (x) lapply(x, mean))


posterior_PIP_cov <- data.frame (pip = unlist(mean_posterior_indicatorcov_by_study), index1 = rep(c(1:36),100), index2 =  rep(c(rep("X1",7),rep("X2",6),rep("X3",3), rep("X4",9), rep("X5",7), rep("X6",4)), 100))

plot_PIP_cov <- ggplot(posterior_PIP_cov, aes(x=index1, y=pip, group = factor(index1), color = factor(index2))) + geom_boxplot()+ylab("Posterior Inclusion Probability")+
                theme_bw()+geom_hline(yintercept=0.5, linetype="dashed", color = "red")+theme(axis.text.y = element_text(size = 10), axis.title.x=element_blank(), axis.title.y=element_text(size = 12),axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(legend.title = element_blank())+ ggtitle("prior inclusion probability ~ Unif(0,1)")


#Median probability model  (PIP >= 0.5)                                                                                               
apply(Result_BSGS_p_unif$Posterior_Inlcusion_Probability_cov,2,mean)

# High frequent model throughout MCMC iterations
ind_posterior_cov <- as.data.frame(Result_BSGS_p_unif$Posterior_Inlcusion_Probability_cov)
rownames(ind_posterior_cov) <- NULL

model_frequency = count(ind_posterior_cov, vars = colnames(ind_posterior_cov))
ind_most_freq_model <- which(model_frequency$freq ==  max(model_frequency$freq))
most_freq_model <- model_frequency[ind_most_freq_model ,]


# Deviance Information Criterion

mean(Result_BSGS_p_unif$DIC, na.rm = T)

#-------------------------------------------------------------------------
# Plots of Figure 4: Posterior Inclusion Probability for levels inclusion
#-------------------------------------------------------------------------

load(file = "Result_BSGS_p_unif.Rdata")
load(file = "Result_BSGS_p_0.5.Rdata")

PIP_levels_0.5 <- split(as.data.frame(Result_BSGS_p_0.5$Posterior_Inlcusion_Probability_within),factor(rep(1:100, each = 3000)))
mean_PIP_levels_0.5 <- lapply(PIP_levels_0.5, function (x) lapply(x, mean))

PIP_levels_unif <- split(as.data.frame(Result_BSGS_p_unif$Posterior_Inlcusion_Probability_within),factor(rep(1:100, each = 3000)))
mean_PIP_levels_unif <- lapply(PIP_levels_unif, function (x) lapply(x, mean))


PIP_levels_all <- data.frame(pip = c(unlist(mean_PIP_levels_0.5), unlist(mean_PIP_levels_unif )), var = rep(1:36,200) , p = rep(c("p = 0.5","p ~ Unif(0,1)"), each = length(unlist(mean_PIP_levels_0.5))))  


#-- X1 --

PIP_levels_all_X1 <- PIP_levels_all[ PIP_levels_all$var %in% 1:7, ]

pip_levels_X1 <- ggplot(PIP_levels_all_X1, aes(var,pip, group = interaction(var,p), fill = p))+geom_boxplot(aes(fill = p, color = p),position=position_dodge(1))+
    ylab("Levels Posterior Inclusion Probability")+geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5)+
     theme(legend.title = element_blank())+theme_bw()+xlab("")+ggtitle(bquote(beta[1] ~ "= (0,1,1,2,2,4,4)"))+ylim(0,1)



                                                                                                
pip_levels_X1 <- pip_levels_X1 +scale_x_continuous(breaks= 1:7,labels=c("level 1", "level 2", "level 3", "level 4", "level 5","level 6", "level 7"))+theme(axis.text.x = element_text(angle = 90, size = 10), axis.text.y = element_text(size = 10))
  

#-- X2 --

PIP_levels_all_X2 <- PIP_levels_all[ PIP_levels_all$var %in% 8:13, ]

pip_levels_X2 <- ggplot(PIP_levels_all_X2, aes(var,pip, group = interaction(var,p)))+ geom_boxplot(aes(fill = p, color = p),position=position_dodge(1))+
  ylab("Levels Posterior Inclusion Probability")+geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5)+
  theme(legend.title = element_blank())+theme_bw()+xlab("")+ggtitle(bquote(beta[2] ~ "= (0,0,0,0,0,0)"))+ylim(0,1)


pip_levels_X2 <- pip_levels_X2 +scale_x_continuous(breaks= 8:13,labels=c("level 1", "level 2", "level 3", "level 4", "level 5","level 6"))+theme(axis.text.x = element_text(angle = 90,size = 10), axis.text.y = element_text(size = 10))

#-- X3 --

PIP_levels_all_X3 <- PIP_levels_all[ PIP_levels_all$var %in% 14:16, ]

pip_levels_X3 <- ggplot(PIP_levels_all_X3, aes(var,pip, group = interaction(var,p)))+ geom_boxplot(aes(fill = p, color = p),position=position_dodge(1))+
  ylab("Levels Posterior Inclusion Probability")+geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5)+
  theme(legend.title = element_blank())+theme_bw()+xlab("")+ggtitle(bquote(beta[3] ~ "= (0,-2,2)"))+ylim(0,1)


pip_levels_X3 <- pip_levels_X3 +scale_x_continuous(breaks= 14:16,labels=c("level 1", "level 2", "level 3"))+theme(axis.text.x = element_text(angle = 90))



#-- X4 --

PIP_levels_all_X4 <- PIP_levels_all[ PIP_levels_all$var %in% 17:25, ]

pip_levels_X4 <- ggplot(PIP_levels_all_X4, aes(var,pip, group = interaction(var,p)))+ geom_boxplot(aes(fill = p, color = p),position=position_dodge(1))+
  ylab("Levels Posterior Inclusion Probability")+geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5)+
  theme(legend.title = element_blank())+theme_bw()+xlab("")+ggtitle(bquote(beta[4] ~ "= (0,0,0,0,0,0,0,0,1)"))+ylim(0,1)


pip_levels_X4 <- pip_levels_X4 +scale_x_continuous(breaks= 17:25,labels=c("level 1", "level 2", "level 3","level 4","level 5","level 6","level 7","level 8","level 9"))+theme(axis.text.x = element_text(angle = 90,size = 10), axis.text.y = element_text(size = 10))

#-- X5 --

PIP_levels_all_X5 <- PIP_levels_all[ PIP_levels_all$var %in% 26:32, ]

pip_levels_X5 <- ggplot(PIP_levels_all_X5, aes(var,pip, group = interaction(var,p)))+ geom_boxplot(aes(fill = p, color = p),position=position_dodge(1))+
  ylab("Levels Posterior Inclusion Probability")+geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5)+
  theme(legend.title = element_blank())+theme_bw()+xlab("")+ggtitle(bquote(beta[5] ~ "= (0,1,1,1,1,-2,-2)"))+ylim(0,1)


pip_levels_X5 <- pip_levels_X5 +scale_x_continuous(breaks= 26:32,labels=c("level 1", "level 2", "level 3","level 4","level 5","level 6","level 7"))+theme(axis.text.x = element_text(angle = 90,size = 10), axis.text.y = element_text(size = 10))

#-- X6 --

PIP_levels_all_X6 <- PIP_levels_all[ PIP_levels_all$var %in% 33:36, ]

pip_levels_X6 <- ggplot(PIP_levels_all_X6, aes(var,pip, group = interaction(var,p)))+ geom_boxplot(aes(fill = p, color = p),position=position_dodge(1))+
  ylab("Levels Posterior Inclusion Probability")+geom_hline(yintercept=0.5, linetype="dashed", color = "red", size = 1.5)+
  theme(legend.title = element_blank())+theme_bw()+xlab("")+ggtitle(bquote(beta[6] ~ "= (0,0,0,0)"))+ylim(0,1)


pip_levels_X6 <- pip_levels_X6 +scale_x_continuous(breaks= 33:36,labels=c("level 1", "level 2", "level 3","level 4"))+theme(axis.text.x = element_text(angle = 90))


