
source('NE_WT_sampler.R')
source('NE_W_moments.R')

NE_STD_sampler = function(n,
                         alpha1,gam1,beta1,
                         alpha2,gam2,beta2,
                         error = 1e-6,init_error = 0.05){
  #n is number of samples wanted
  #theta is probability negative, decided by E[X_1]/(E[X_1]- E[X_2]) 
  #param_1 => parameters for >0
  #param_2 => parameters for <0
  #error and init error are numerical scheme parameters
  
  
  
  #Calculate theta
  mean_1 = NE_W_moment(1,c(0,Inf),
                       alpha = alpha1, gam = gam1,
                       beta = beta1)[[1]]
  mean_2 = -NE_W_moment(1,c(0,Inf),
                        alpha = alpha2, gam = gam2,
                        beta = beta2)[[1]]
  theta = mean_1/(mean_1-mean_2)
  print(theta)
  return(NE_WT_sampler(n,theta,
                       alpha1,gam1,beta1,
                       alpha2,gam2,beta2,)
         )
  

}


#TESTER
# NE_STD_var1_sampler = function(n,
#                           alpha1,gam1,beta1,
#                           alpha2,gam2,beta2,
#                           error = 1e-6,init_error = 0.05){
#   #n is number of samples wanted
#   #theta is probability negative, decided by E[X_1]/(E[X_1]- E[X_2]) 
#   #param_1 => parameters for >0
#   #param_2 => parameters for <0
#   #error and init error are numerical scheme parameters
#   
#   
#   
#   #Calculate theta
#   mean_1 = NE_W_moment(1,c(0,Inf),
#                        alpha = alpha1, gam = gam1,
#                        beta = beta1)[[1]]
#   mean_2 = -NE_W_moment(1,c(0,Inf),
#                         alpha = alpha2, gam = gam2,
#                         beta = beta2)[[1]]
#   theta = mean_1/(mean_1-mean_2)
#   
#   
#   var_ne = moment_func(2, NE_WT_PDF, domain = c(-Inf,Inf),
#                          theta = theta,
#                          alpha1 = alpha1,gam1 = gam1,beta1 = beta1,
#                          alpha2 = alpha2,gam2 = gam2,beta2 = beta2)[[1]]
#   samples = NE_WT_sampler(n,theta,
#                           alpha1,gam1,beta1,
#                           alpha2,gam2,beta2,)
#   samples[1,] = samples[1,]/sqrt(var_ne)
#   return(samples)
#   
#   
# }
# 
