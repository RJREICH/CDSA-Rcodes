#------------------------------------------------------------
# Function for Bayesian effect fusion model based clustering 
#------------------------------------------------------------

# Data_simulated_nstudies: Simulated data obtained from Simulated_data R script
# nb_studies: nb of silulated studies
# nb_iter: number of MCMC iterations
#nu: hyperparameter to set by the user
#e: hyperparameter to set by the user

BEF <- function(Data_simulated_nstudies,nb_studies,nb_iter,nu,e){
  
  #load packages
  library(doParallel)
  library(foreach)
  
  cores <- detectCores()
  # create clusters
  cl <- makeCluster(cores-3)
  registerDoParallel(cl)
  
  #Import JAGS BSGL-SS
  
  model.jags <- paste(getwd(),"/BEF.txt", sep="")
  
  result = foreach (i=1:nb_studies,.packages = c("R2jags"),.combine = rbind)%dopar% {
    
    set.seed(123)
    
    # The ith simulated data set   
    Data_sim <- Data_simulated_nstudies[[i]]
    
    # Hyperparameter specifications (Malsiner et al., 2017)
    X = as.matrix(Data_sim[,-1])
    lm.model <- lm(Data_sim$y~ X-1)
    
    #Regression coeffcients estimations
    Beta1 <- as.numeric(coefficients(lm.model)[2:8])
    Beta2 <- as.numeric(coefficients(lm.model)[9:14])
    Beta3 <- as.numeric(coefficients(lm.model)[15:17])
    Beta4 <- as.numeric(coefficients(lm.model)[18:26])
    Beta5 <- as.numeric(coefficients(lm.model)[27:33])
    Beta6 <- as.numeric(coefficients(lm.model)[34:37])
    
    #Precision estimations
    tau <- 1/summary(lm.model)$sigma^2
    
    # choice of Hyperparameters  
    
    m_B1=mean(Beta1)
    m_B2=mean(Beta2)
    m_B3=mean(Beta3)
    m_B4=mean(Beta4)
    m_B5=mean(Beta5)
    m_B6=mean(Beta6)
    
    inv_M0_1 <- 1/(min(Beta1)-max(Beta1))^2;
    inv_M0_2 <- 1/(min(Beta2)-max(Beta2))^2;
    inv_M0_3 <- 1/(min(Beta3)-max(Beta3))^2;
    inv_M0_4 <- 1/(min(Beta4)-max(Beta4))^2;
    inv_M0_5 <- 1/(min(Beta5)-max(Beta5))^2;
    inv_M0_6 <- 1/(min(Beta6)-max(Beta6))^2;
    
    
    V1 <- (1/(7-1))*sum((Beta1-m_B1)^2)
    V2 <- (1/(6-1))*sum((Beta2-m_B2)^2)
    V3 <- (1/(3-1))*sum((Beta3-m_B3)^2)
    V4 <- (1/(9-1))*sum((Beta4-m_B4)^2)
    V5 <- (1/(7-1))*sum((Beta5-m_B5)^2)
    V6 <- (1/(4-1))*sum((Beta6-m_B6)^2)
    V = c(V1,V2,V3,V4,V5,V6)
    
    # JAGS Data
    
    data = list( e_1 = rep(e,8), e_2 = rep(e,7),e_3 = rep(e,4),
                 e_4 = rep(e,10), e_5 = rep(e,8), e_6 = rep(e,5),
                 N = nrow(X), y=Data_sim$y, X =Data_sim[,-1],level = c(7,6,3,9,7,4),I_1 = diag(7),
                 I_2 = diag(6) , I_3 = diag(3) , I_4= diag(9), I_5 = diag(7) ,I_6 = diag(4),
                 V= as.numeric(V), nu = nu,
                 inv_M0_1 = inv_M0_1,inv_M0_2 = inv_M0_2,inv_M0_3 = inv_M0_3,inv_M0_4 = inv_M0_4,inv_M0_5 = inv_M0_5,inv_M0_6 = inv_M0_6,
                 mean_B1 = rep(m_B1,7),mean_B2 = rep(m_B2,6),mean_B3 = rep(m_B3,3),mean_B4 = rep(m_B4,9),mean_B5 = rep(m_B5,7),
                 mean_B6 = rep(m_B6,4))
    
    # JAGS initialization
    
    inits <-  list(list(BB1 = Beta1 , BB2 = Beta2, BB3 = Beta3, BB4 = Beta4, BB5 = Beta5, BB6 = Beta6, tau = tau ),
                   list(BB1 = runif(1,0.9,1)*Beta1 , BB2 = runif(1,0.9,1)*Beta2, BB3 = runif(1,0.9,1)*Beta3, BB4 = runif(1,0.9,1)*Beta4, BB5 = runif(1,0.9,1)*Beta5, BB6 = runif(1,0.9,1)*Beta6, tau = runif(1,0.9,1)*tau ),
                   list(BB1 = runif(1,0.9,1)*Beta1 , BB2 = runif(1,0.9,1)*Beta2, BB3 = runif(1,0.9,1)*Beta3, BB4 = runif(1,0.9,1)*Beta4, BB5 = runif(1,0.9,1)*Beta5, BB6 = runif(1,0.9,1)*Beta6, tau = runif(1,0.9,1)*tau ))
    
    
    # JAGS traced parameters
    
    params <- c("B", "tau","lambda1","lambda2","lambda3","lambda4","lambda5","lambda6","inv_psi","z1","z2","z3","z4","z5","z6")
    
    
    # Run JAGS model
    
    
    jags_fusion_model <- jags(data = data,parameters.to.save = params,n.iter = nb_iter
                             ,model.file =  model.jags,n.chains = 3, inits =inits)
    
    # loop output
    cbind(jags_fusion_model$BUGSoutput$sims.list$B, jags_fusion_model$BUGSoutput$sims.list$z1, jags_fusion_model$BUGSoutput$sims.list$z2,
          jags_fusion_model$BUGSoutput$sims.list$z3, jags_fusion_model$BUGSoutput$sims.list$z4,
          jags_fusion_model$BUGSoutput$sims.list$z5, jags_fusion_model$BUGSoutput$sims.list$z6, c(jags_fusion_model$BUGSoutput$summary[,8], rep(NA, length(jags_fusion_model$BUGSoutput$sims.list$tau)-length(jags_fusion_model$BUGSoutput$summary[,8]))), c(jags_fusion_model$BUGSoutput$DIC, rep(NA,length(jags_fusion_model$BUGSoutput$sims.list$tau)-1)))
    
  }
  # stop clusters
  stopCluster(cl)
  
  list(Posterior_Effect_Regression = result[,1:37]  , 
       Posterior_Inlcusion_Probability_z1 = result[,38:45],
       Posterior_Inlcusion_Probability_z2 = result[,46:52],
       Posterior_Inlcusion_Probability_z3 = result[,53:56],
       Posterior_Inlcusion_Probability_z4 = result[,57:66],
       Posterior_Inlcusion_Probability_z5 = result[,67:74],
       Posterior_Inlcusion_Probability_z6 = result[,75:79],
       Rhat = result[,80], 
       DIC = result[,81])}
  
#----------------------------  
# Application of BEF function
#----------------------------

# Functions inputs
load(file = "Data_sim_100obs_100studies.Rdata")

#-------------------
# e = 0.01, nu = 100
#-------------------

T1.BEF_e_0.01_nu_100 <- Sys.time()
Result_BEF_e_0.01_nu_100 = BEF(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 100000, e = 0.01, nu = 100)
T2.BEF_e_0.01_nu_100 <- Sys.time()

Time.BEF_e_0.01_nu_100 <- T2.BEF_e_0.01_nu_100-T1.BEF_e_0.01_nu_100
save(Time.BEF_e_0.01_nu_100, file = "Time.BEF_e_0.01_nu_100.Rdata")  
save(Result_BEF_e_0.01_nu_100, file = "Result_BEF_e_0.01_nu_100.Rdata")


#-------------------
# e = 0.1, nu = 100
#-------------------

T1.BEF_e_0.1_nu_100 <- Sys.time()
Result_BEF_e_0.1_nu_100 = BEF(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 100000, e = 0.1, nu = 100)
T2.BEF_e_0.1_nu_100 <- Sys.time()

Time.BEF_e_0.1_nu_100 <- T2.BEF_e_0.1_nu_100-T1.BEF_e_0.1_nu_100
save(Time.BEF_e_0.1_nu_100, file = "Time.BEF_e_0.1_nu_100.Rdata")
save(Result_BEF_e_0.1_nu_100, file = "Result_BEF_e_0.1_nu_100.Rdata")

#-------------------
# e = 0.6, nu = 100
#-------------------

T1.BEF_e_0.6_nu_100 <- Sys.time()
Result_BEF_e_0.6_nu_100 = BEF(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 100000, e = 0.6, nu = 100)
T2.BEF_e_0.6_nu_100 <- Sys.time()

Time.BEF_e_0.6_nu_100 <- T2.BEF_e_0.6_nu_100-T1.BEF_e_0.6_nu_100
save(Time.BEF_e_0.6_nu_100, file = "Time.BEF_e_0.6_nu_100.Rdata")
save(Result_BEF_e_0.6_nu_100, file = "Result_BEF_e_0.6_nu_100.Rdata")
#-------------------
# e = 0.01, nu = 10
#-------------------

T1.BEF_e_0.01_nu_10 <- Sys.time()
Result_BEF_e_0.01_nu_10 = BEF(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 100000, e = 0.01, nu = 10)
T2.BEF_e_0.01_nu_10 <- Sys.time()

Time.BEF_e_0.01_nu_10 <- T2.BEF_e_0.01_nu_10-T1.BEF_e_0.01_nu_10
save(Time.BEF_e_0.01_nu_10, file = "Time.BEF_e_0.01_nu_10.Rdata")
save(Result_BEF_e_0.01_nu_10, file = "Result_BEF_e_0.01_nu_10.Rdata")

#---------------------
# e = 0.01, nu = 1000
#---------------------


T1.BEF_e_0.01_nu_1000 <- Sys.time()
Result_BEF_e_0.01_nu_1000 = BEF(Data_simulated_nstudies = Data_sim_100obs_100studies,nb_studies = 100,nb_iter = 100000, e = 0.01, nu = 1000)
T2.BEF_e_0.01_nu_1000 <- Sys.time()

Time.BEF_e_0.01_nu_1000 <- T2.BEF_e_0.01_nu_1000-T1.BEF_e_0.01_nu_1000
save(Time.BEF_e_0.01_nu_1000, file = "Time.BEF_e_0.01_nu_1000.Rdata")
save(Result_BEF_e_0.01_nu_1000, file = "Result_BEF_e_0.01_nu_1000.Rdata")
 
#----------------------------------------------------------------------
# Plots of Figure 5: Bayesian Effect Fusion with e = 0.01 and nu = 100
#----------------------------------------------------------------------

# load packages
library(reshape2)
library(ggplot2)
library(ggcorrplot)

# The posterior levels fusion
ind_var_1 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z1
# > ind_var_1[1,]
# [1] 8 5 5 5 5 6 6 6 which mean that level 2, level 3, level 4, level 5 are fused together and form one new level
# in addition, level 6, level 7, level 8 belongs to the same cluster
ind_var_2 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z2
ind_var_3 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z3
ind_var_4 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z4
ind_var_5 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z5
ind_var_6 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z6

# counting how much levels are fused together among all iterations


M1 <- array(0, dim = c(8,8,nrow(ind_var_1)))  # 8 is the number of X1 levels+1
M2 <- array(0, dim = c(7,7, nrow(ind_var_2)))
M3 <- array(0, dim = c(4,4,nrow(ind_var_3)))
M4 <- array(0,dim = c(10,10, nrow(ind_var_4)))
M5 <- array(0,dim = c(8,8, nrow(ind_var_5)))
M6 <- array(0,dim = c(5,5, nrow(ind_var_6)))


for (i in 1:nrow(ind_var_1)){
  
  m1 = ind_var_1[i,]
  m2 = ind_var_2[i,]
  m3 = ind_var_3[i,]
  m4 = ind_var_4[i,]
  m5 = ind_var_5[i,]
  m6 = ind_var_6[i,]
  
  M1[,,i]= 1*outer(m1,m1, "==")
  M2[,,i]= 1*outer(m2,m2, "==")
  M3[,,i]= 1*outer(m3,m3, "==")
  M4[,,i]= 1*outer(m4,m4, "==")
  M5[,,i]= 1*outer(m5,m5, "==")
  M6[,,i]= 1*outer(m6,m6, "==")
}

# Computing levels posterior mean fusion for each of categorical covariate

#-- X1--
mean_F1 <-  apply(M1,1:2,mean)
row.names(mean_F1) <- colnames(mean_F1) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7","level 0")


#-- X2--

mean_F2 <-  apply(M2,1:2,mean)
row.names(mean_F2) <- colnames(mean_F2) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 0")

#-- X3--

mean_F3 <-  apply(M3,1:2,mean)
row.names(mean_F3) <- colnames(mean_F3) <- c("level 1","level 2","level 3", "level 0")

#-- X4--

mean_F4 <-  apply(M4,1:2,mean)
row.names(mean_F4) <- colnames(mean_F4) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7","level 8","level 9", "level 0")

#-- X5--

mean_F5 <-  apply(M5,1:2,mean)
row.names(mean_F5) <- colnames(mean_F5) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7", "level 0")


#-- X6--

mean_F6 <-  apply(M6,1:2,mean)
row.names(mean_F6) <- colnames(mean_F6) <- c("level 1","level 2","level 4","level 5", "level 0")


# Matrix plots for each of the categorical covariate

plot_F1 <- ggcorrplot(mean_F1, lab = T, title = bquote(beta[1] ~ "= (0,1,1,2,2,4,4)"))+theme(legend.position='none')
plot_F2 <- ggcorrplot(mean_F2, lab = T,title = bquote(beta[2] ~ "= (0,0,0,0,0,0)"))+theme(legend.position='none')
plot_F3 <- ggcorrplot(mean_F3, lab = T,title = bquote(beta[3] ~ "= (0,-2,2)"))+theme(legend.position='none')
plot_F4 <- ggcorrplot(mean_F4, lab = T,title = bquote(beta[4] ~ "= (0,0,0,0,0,0,0,0,1)"))+theme(legend.position='none')
plot_F5 <- ggcorrplot(mean_F5, lab = T,title = bquote(beta[5] ~ "= (0,1,1,1,1,-2,-2)"))+theme(legend.position='none')
plot_F6 <- ggcorrplot(mean_F6, lab = T,title = bquote(beta[6] ~ "= (0,0,0,0)"))+theme(legend.position='none')

#---------------------------------------------------------------------------------
# Code for Figure 6: Median of the posterior regression effects for all covariates
# based on Bayesian Effect Fusion with e = 0.01 and nu = 100
#---------------------------------------------------------------------------------


B_posterior_estimation <- Result_BEF_e_0.01_nu_100$Posterior_Effect_Regression
B_median_estimation <- apply(B_posterior_estimation,2,median)



Median_all <- data.frame (m = B_median_estimation, index2 = c("intercept",rep("X1",7),rep("X2",6),rep("X3",3), rep("X4",9), rep("X5",7), rep("X6",4)))


plot_median_all <- ggplot(Median_all, aes(x=1:37, y=m, group = factor(index2), color = factor(index2))) + geom_point(size = 2.5)+ xlab("")+ylab("Posterior Median Regression Effects")+
                   geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.7)+theme_bw()+theme(legend.title = element_blank())+theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 10), axis.ticks = element_blank())

plot_median_all

#------------------------------------------------------------
# DIC comparaison model for e = 0.01, 0.1 and 0.6 in Table 4
#------------------------------------------------------------

load(file = "Result_BEF_e_0.01_nu_100.Rdata")
load(file = "Result_BEF_e_0.1_nu_100.Rdata")
load(file = "Result_BEF_e_0.6_nu_100.Rdata")

DIC_e_0.01_nu_100 <- mean(Result_BEF_e_0.01_nu_100$DIC, na.rm = T)
DIC_e_0.1_nu_100 <- mean(Result_BEF_e_0.1_nu_100$DIC, na.rm = T)
DIC_e_0.6_nu_100 <- mean(Result_BEF_e_0.6_nu_100$DIC, na.rm = T)

#----------------------------------------
# Matrix plots in Figure 7 for X2 and X4
#----------------------------------------


load(file = "Result_BEF_e_0.01_nu_100.Rdata")
load(file = "Result_BEF_e_0.1_nu_100.Rdata")
load(file = "Result_BEF_e_0.6_nu_100.Rdata")


# The posterior levels fusion

# e = 0.01, nu = 100
ind_var_X2_e_0.01 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z2
ind_var_X4_e_0.01 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z4

# e = 0.1, nu = 100
ind_var_X2_e_0.1 <- Result_BEF_e_0.1_nu_100$Posterior_Inlcusion_Probability_z2
ind_var_X4_e_0.1 <- Result_BEF_e_0.1_nu_100$Posterior_Inlcusion_Probability_z4

# e = 0.6, nu = 100
ind_var_X2_e_0.6 <- Result_BEF_e_0.6_nu_100$Posterior_Inlcusion_Probability_z2
ind_var_X4_e_0.6 <- Result_BEF_e_0.6_nu_100$Posterior_Inlcusion_Probability_z4


# counting how much levels are fused together among all iteration

M1 <- array(0,dim = c(7,7,nrow(ind_var_X2_e_0.01)))  # 6 is the number of X2 levels
M2 <- array(0,dim = c(10,10, nrow(ind_var_X4_e_0.01))) # 9 is the number of X4 levels
M3 <- array(0,dim = c(7,7,nrow(ind_var_X2_e_0.1)))
M4 <- array(0,dim = c(10,10, nrow(ind_var_X4_e_0.1)))
M5 <- array(0,dim = c(7,7, nrow(ind_var_X2_e_0.6)))
M6 <- array(0,dim = c(10,10, nrow(ind_var_X4_e_0.6)))

for (i in 1:nrow(ind_var_X2_e_0.01)){
  
  m1 = ind_var_X2_e_0.01[i,]
  m2 = ind_var_X4_e_0.01[i,]
  m3 = ind_var_X2_e_0.1[i,]
  m4 = ind_var_X4_e_0.1[i,]
  m5 = ind_var_X2_e_0.6[i,]
  m6 = ind_var_X4_e_0.6[i,]
  
  M1[,,i]= 1*outer(m1,m1, "==")
  M2[,,i]= 1*outer(m2,m2, "==")
  M3[,,i]= 1*outer(m3,m3, "==")
  M4[,,i]= 1*outer(m4,m4, "==")
  M5[,,i]= 1*outer(m5,m5, "==")
  M6[,,i]= 1*outer(m6,m6, "==")
}

# Computing levels posterior mean fusion for each of categorical covariate

# panels 1 and 4 of Figure 7

mean_X2_e_0.01_nu_100 <-  apply(M1,1:2,mean)
row.names(mean_X2_e_0.01_nu_100) <- colnames(mean_X2_e_0.01_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6", "level 0")


mean_X4_e_0.01_nu_100 <-  apply(M2,1:2,mean)
row.names(mean_X4_e_0.01_nu_100) <- colnames(mean_X4_e_0.01_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7","level 8","level 9", "level 0")

plot_X2_e_0.01_nu_100 <- ggcorrplot(mean_X2_e_0.01_nu_100, lab = T,title = bquote(beta[2] ~ "= (0,0,0,0,0,0)"))+theme(legend.position='none')
plot_X4_e_0.01_nu_100 <- ggcorrplot(mean_X4_e_0.01_nu_100, lab = T,title = bquote(beta[4] ~ "= (0,0,0,0,0,0,0,0,1)"))+theme(legend.position='none')


# panels 2 and 5 of Figure 7

mean_X2_e_0.1_nu_100 <-  apply(M3,1:2,mean)
row.names(mean_X2_e_0.1_nu_100) <- colnames(mean_X2_e_0.1_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6", "level 0")


mean_X4_e_0.1_nu_100 <-  apply(M4,1:2,mean)
row.names(mean_X4_e_0.1_nu_100) <- colnames(mean_X4_e_0.1_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7","level 8","level 9", "level 0")

plot_X2_e_0.1_nu_100 <- ggcorrplot(mean_X2_e_0.1_nu_100, lab = T,title = bquote(beta[2] ~ "= (0,0,0,0,0,0)"))+theme(legend.position='none')
plot_X4_e_0.1_nu_100 <- ggcorrplot(mean_X4_e_0.1_nu_100, lab = T,title = bquote(beta[4] ~ "= (0,0,0,0,0,0,0,0,1)"))+theme(legend.position='none')

# panels 3 and 6 of Figure 7

mean_X2_e_0.6_nu_100 <-  apply(M5,1:2,mean)
row.names(mean_X2_e_0.6_nu_100) <- colnames(mean_X2_e_0.6_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6", "level 0")


mean_X4_e_0.6_nu_100 <-  apply(M6,1:2,mean)
row.names(mean_X4_e_0.6_nu_100) <- colnames(mean_X4_e_0.6_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7","level 8","level 9", "level 0")

plot_X2_e_0.6_nu_100 <- ggcorrplot(mean_X2_e_0.6_nu_100, lab = T,title = bquote(beta[2] ~ "= (0,0,0,0,0,0)"))+theme(legend.position='none')
plot_X4_e_0.6_nu_100 <- ggcorrplot(mean_X4_e_0.6_nu_100, lab = T,title = bquote(beta[4] ~ "= (0,0,0,0,0,0,0,0,1)"))+theme(legend.position='none')


#-----------------------------------------------------------
# DIC comparaison model for nu = 10, 100 and 1000 in Table 5
#-----------------------------------------------------------

load(file = "Result_BEF_e_0.01_nu_10.Rdata")
load(file = "Result_BEF_e_0.01_nu_100.Rdata")
load(file = "Result_BEF_e_0.01_nu_1000.Rdata")

DIC_e_0.01_nu_10 <- mean(Result_BEF_e_0.01_nu_10$DIC, na.rm = T)
DIC_e_0.01_nu_100 <- mean(Result_BEF_e_0.01_nu_100$DIC, na.rm = T)
DIC_e_0.01_nu_1000 <- mean(Result_BEF_e_0.01_nu_1000$DIC, na.rm = T)

#---------------------------------------
# Matrix plots in Figure 8 for X1 and X5
#---------------------------------------


load(file = "Result_BEF_e_0.01_nu_10.Rdata")
load(file = "Result_BEF_e_0.01_nu_100.Rdata")
load(file = "Result_BEF_e_0.01_nu_1000.Rdata")


# The posterior levels fusion

# e = 0.01, nu = 10
ind_var_X1_nu_10 <- Result_BEF_e_0.01_nu_10$Posterior_Inlcusion_Probability_z1[,-8]
ind_var_X5_nu_10 <- Result_BEF_e_0.01_nu_10$Posterior_Inlcusion_Probability_z5[,-8]

# e = 0.01, nu = 100
ind_var_X1_nu_100 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z1[,-8]
ind_var_X5_nu_100 <- Result_BEF_e_0.01_nu_100$Posterior_Inlcusion_Probability_z5[,-8]

# e = 0.01, nu = 1000
ind_var_X1_nu_1000 <- Result_BEF_e_0.01_nu_1000$Posterior_Inlcusion_Probability_z1[,-8]
ind_var_X5_nu_1000 <- Result_BEF_e_0.01_nu_1000$Posterior_Inlcusion_Probability_z5[,-8]


# counting how much levels are fused together among all iterations

M1 <- array(0,dim = c(7,7,nrow(ind_var_X1_nu_10)))  # 7 is the number of X1 levels
M2 <- array(0,dim = c(7,7,nrow(ind_var_X5_nu_10))) # 7 is the number of X5 levels
M3 <- array(0,dim = c(7,7,nrow(ind_var_X1_nu_100)))
M4 <- array(0,dim = c(7,7,nrow(ind_var_X5_nu_100)))
M5 <- array(0,dim = c(7,7,nrow(ind_var_X1_nu_1000)))
M6 <- array(0,dim = c(7,7,nrow(ind_var_X5_nu_1000)))

for (i in 1:nrow(ind_var_X1_nu_10)){
  
  m1 = ind_var_X1_nu_10[i,]
  m2 = ind_var_X5_nu_10[i,]
  m3 = ind_var_X1_nu_100[i,]
  m4 = ind_var_X5_nu_100[i,]
  m5 = ind_var_X1_nu_1000[i,]
  m6 = ind_var_X5_nu_1000[i,]
  
  M1[,,i]= 1*outer(m1,m1, "==")
  M2[,,i]= 1*outer(m2,m2, "==")
  M3[,,i]= 1*outer(m3,m3, "==")
  M4[,,i]= 1*outer(m4,m4, "==")
  M5[,,i]= 1*outer(m5,m5, "==")
  M6[,,i]= 1*outer(m6,m6, "==")
}

# Computing levels posterior mean fusion foreach of categorical covariate

# panels 1 and 5 of Figure 8

mean_X1_e_0.01_nu_10 <-  apply(M1,1:2,mean)
row.names(mean_X1_e_0.01_nu_10) <- colnames(mean_X1_e_0.01_nu_10) <- c("level 1","level 2","level 3","level 4","level 5","level 6", "level 7")


mean_X5_e_0.01_nu_10 <-  apply(M2,1:2,mean)
row.names(mean_X5_e_0.01_nu_10) <- colnames(mean_X5_e_0.01_nu_10) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7")

plot_X1_e_0.01_nu_10 <- ggcorrplot(mean_X1_e_0.01_nu_10, lab = T,title = bquote(beta[1] ~ "= (0,1,1,2,2,4,4)"))+theme(legend.position='none')
plot_X5_e_0.01_nu_10<- ggcorrplot(mean_X5_e_0.01_nu_10, lab = T,title = bquote(beta[5] ~ "= (0,1,1,1,1,-2,-2)"))+theme(legend.position='none')


# panels 2 and 5 of Figure 8

mean_X1_e_0.01_nu_100 <-  apply(M3,1:2,mean)
row.names(mean_X1_e_0.01_nu_100) <- colnames(mean_X1_e_0.01_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6", "level 7")


mean_X5_e_0.01_nu_100 <-  apply(M4,1:2,mean)
row.names(mean_X5_e_0.01_nu_100) <- colnames(mean_X5_e_0.01_nu_100) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7")

plot_X1_e_0.01_nu_100 <- ggcorrplot(mean_X1_e_0.01_nu_100, lab = T,title = bquote(beta[1] ~ "= (0,1,1,2,2,4,4)"))+theme(legend.position='none')
plot_X5_e_0.01_nu_100 <- ggcorrplot(mean_X5_e_0.01_nu_100, lab = T,title = bquote(beta[5] ~ "= (0,1,1,1,1,-2,-2)"))+theme(legend.position='none')


# panels 3 and 6 of Figure 8

mean_X1_e_0.01_nu_1000 <-  apply(M5,1:2,mean)
row.names(mean_X1_e_0.01_nu_1000) <- colnames(mean_X1_e_0.01_nu_1000) <- c("level 1","level 2","level 3","level 4","level 5","level 6", "level 7")


mean_X5_e_0.01_nu_1000 <-  apply(M6,1:2,mean)
row.names(mean_X5_e_0.01_nu_1000) <- colnames(mean_X5_e_0.01_nu_1000) <- c("level 1","level 2","level 3","level 4","level 5","level 6","level 7")

plot_X1_e_0.01_nu_1000 <- ggcorrplot(mean_X1_e_0.01_nu_1000, lab = T,title = bquote(beta[1] ~ "= (0,1,1,2,2,4,4)"))+theme(legend.position='none')
plot_X5_e_0.01_nu_1000 <- ggcorrplot(mean_X5_e_0.01_nu_1000, lab = T,title = bquote(beta[5] ~ "= (0,1,1,1,1,-2,-2)"))+theme(legend.position='none')


