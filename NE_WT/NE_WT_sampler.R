
source('NE_W_sampler.R')

NE_WT_sampler = function(n,theta,
                         alpha1,gam1,beta1,
                         alpha2,gam2,beta2,
                         error = 1e-6,init_error = 0.05){
  #n is number of samples wanted
  #theta is probability negative 
  #param_1 => parameters for >0
  #param_2 => parameters for <0
  #error and init error are numerical scheme parameters
  
  #sample on latent Z which decides how many negative and positive
  Z = rbinom(n,1,theta) + 1
    #Z_i = 2 => x_i is negative;  Z_i = 1 => x_i is positive
  
  n_pos = length(which(Z == 1))
  pos_idx = which(Z==1)
  n_neg = n-n_pos
  neg_idx = which(Z == 2)
  
  x_pos = NE_W_sampler(n_pos,alpha = alpha1, gam = gam1,
                       beta = beta1)
  x_pos[,2] = x_pos[,2]*(1-theta) + theta
  x_neg = -NE_W_sampler(n_neg,alpha = alpha2, gam = gam2,
                       beta = beta2)
  x_neg[,2] = x_neg[,2]*theta + theta
  
  x = matrix(rep.int(0, 2*n),ncol = 2)

  x[pos_idx,] = x_pos
  x[neg_idx,] = x_neg
  return(x)
  
}

