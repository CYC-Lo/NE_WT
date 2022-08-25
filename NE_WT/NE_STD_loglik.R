# Log Likelihood equation for NE_WT

#Due to optim requirements, it takes all parameters(vector) as first argument
source('NE_W_loglik.R')
source('NE_W_moments.R')



NE_STD_logLik = function(params,x){
  # params is in order alpha1,gam1,beta1,alpha2,gam2,beta2
  
  x_1 = x[x >= 0]
  x_2 = -(x[x < 0])
  N_1 = length(x_1)
  N_2 = length(x_2)
  params_1 = params[1:3]
  params_2 = params[4:6]
  loglik_1 = NE_W_logLik(params_1,x_1)
  loglik_2 = NE_W_logLik(params_2,x_2)

  
  mean_1 = NE_W_moment(1,c(0,Inf),
                       alpha = params_1[1], gam = params_1[2],
                       beta = params_1[3])[[1]]
  mean_2 = -NE_W_moment(1,c(0,Inf),
                       alpha = params_2[1], gam = params_2[2],
                       beta = params_2[3])[[1]]
  theta = mean_1/(mean_1-mean_2)
  loglik = N_1*log(1-theta) + N_2*log(theta) + loglik_1+loglik_2
  return(loglik)
}


#Obsolete theta section
# NE_STD_theta_partial_logLik = function(theta, N_1, N_2){
#   # params is significantly different to other partials
#   theta = params[1]
#   L = N_2/theta - N_1/(1-theta)
#   return(L)
# }



NE_STD_partial_loglik = function(params,x){
  # params is in order alpha1,gam1,beta1,alpha2,gam2,beta2
  #theta is not included since it's known from E[X_1]/(E[X_1]-E[X_2]
  
  x_1 = x[x >= 0]
  x_2 = -(x[x < 0])
  # x_1 is postive real line data
  # x_2 is  absolute value of negatice real line data
  
  N_1 = length(x_1)
  N_2 = length(x_2)
  
  params_1 = params[1:3]
  params_2 = params[4:6]
  
  # Obsolete theta section
  # mean_1 = NE_W_moment(1,c(0,Inf),
  #                      alpha = params_1[1], gam = params_1[2],
  #                      beta = params_1[3])[[1]]
  # mean_2 = -NE_W_moment(1,c(0,Inf),
  #                       alpha = params_2[1], gam = params_2[2],
  #                       beta = params_2[3])[[1]]
  # theta = mean_1/(mean_1-mean_2)
  # 

  
  partial_1 = NE_W_partial_loglik(params_1,x_1)
  partial_2= NE_W_partial_loglik(params_2,x_2)
  return(c( #NE_STD_theta_partial_logLik(theta,N_1,N_2),
           partial_1,
           partial_2)
         )
  
}

NE_STD_MLE = function(data, init_params = c(2,2,2,2,2,2),...){
  Results = optim(par = init_params, fn = NE_STD_logLik,
                  gr = NE_STD_partial_loglik, x = data,
                  control = list(fnscale = -1)
  ) 
  
}



