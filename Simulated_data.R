#----------------------------------------
# Function for simulated data generation
#----------------------------------------

# nb_obs: number of observations
# nb_studies: number of simulated studies
# nb_categorical_cov: number of categorical covariates
# vec_nb_levels: number of levels with each categorical covariate
# B.vector: regression coefficients
# sigma_y: standard deviation of the linear model 


Data_simulation <- function(nb_obs,nb_studies,nb_categorical_cov, vec_nb_levels,B.vector,sigma_y){
  set.seed(123)
  
  # Create empty list to store data frames
  Data_sim <-  vector("list", length = nb_studies)
  
  for (i in 1:nb_studies){

  # Generate randomly nb_categorical_cov categorical covariates with vec_nb_levels respectively
  
  # X is defiend as the design matrix (each column represents a covariate)
     X <- data.frame(matrix(NA,ncol=nb_categorical_cov,nrow=nb_obs))
  
     for (j in 1:nb_categorical_cov){
       X[,j] <- as.factor(sample(1:vec_nb_levels[j],nb_obs, replace = T))
     }
  
  # Expanding factors to a set of dummy varaibles using treatment contrast
  
  X.design <- model.matrix(~.,data = as.data.frame(X))
  
  # Generating the response y from a linear model 
  
  y <- rnorm(nb_obs,X.design%*%B.vector,sigma_y)
  
  # Stock the generated data set as a list element of Data_sim
  Data_sim[[i]] <- data.frame(y = y, X.design)
  }
  
  Data_sim
  
  }

#----------------------------------------------------------------
# Generate artficial data corresponding to the following features
#----------------------------------------------------------------

nb_obs = 100
nb_studies = 100
nb_categorical_cov = 6
vec_nb_levels = c(8,7,4,10,8,5)
B.vector = c(1,0,1,1,2,2,4,4,0,0,0,0,0,0,0,-2,2,0,0,0,0,0,0,0,0,1,0,1,1,1,1,-2,-2,0,0,0,0)
sigma_y = 0.5

Data_sim_100obs_100studies = Data_simulation(nb_obs = nb_obs,nb_studies = nb_studies,nb_categorical_cov = nb_categorical_cov, vec_nb_levels = vec_nb_levels,B.vector = B.vector,sigma_y = sigma_y)

#-------------
# save results
#-------------

save(Data_sim_100obs_100studies, file = "Data_sim_100obs_100studies.Rdata")  




