#Alternative Log Likelihood function for NE-WGARCH
#Same actual function
#Different constrainting parameter function

source('../NE_WT/NE_STD_loglik.R',chdir = TRUE)


NE_WGARCH_loglik_alt = function(params,r,r0){
  # params is in order alpha1,gam1,beta1,alpha2,gam2,beta2,
  #                    pre_w,pre_a,pre_b,sigma0
  
  # Basically NE_WT params then GARCH params
  
  # r is data, r or epsilon from 1:T (not the x from notation)
  # r0 is epsilon 0, the data point before first day of interest
  
  # NE_WT params 
  #Need to transform to ensure all > 0
  #   as well as alpha*beta -1 \geq 0
  params_1 = abs(params[1:3])
  params_1[1] = params_1[1] + params_1[3]^(-1)
  params_2 = abs(params[4:6])
  params_2[1] = params_2[1] + params_2[3]^(-1)
  # GARCH params
  #Need to transform to fit constraint
  w = abs(params[7])
  a = 0.5*(params[8]/sqrt(1+params[8]^2) +1 )
  b = (1-a)*0.5*(params[9]/sqrt(1+params[9]^2) +1 )
  sigma0 = abs(params[10])
  
  
  T_len = length(r)
  
  
  # step 1
  sigma_sq = c()
  sigma_sq[1] = w + a*r0^2 + b*sigma0^2
  for (i in 2:T_len) {
    sigma_sq[i] = w + a*r[i-1]^2 + b*sigma_sq[i-1]
  }
  
  
  # step 2
  x = r/sqrt(sigma_sq)
  # return(list(sigma_sq = sigma_sq, x = x))
  # step 3
  L_NE = NE_STD_logLik(params = c(params_1,params_2), x = x)
  L_sig = sum(0.5*log(sigma_sq)) # log(sqrt) = 0.5*log
  return(L_NE - L_sig)
  
}

NE_WGARCH_MLE_alt = function(data, init_params = c(1,1,1,1,1,1,
                                               1,1,2,0.5), maxit = 1000,...){
  #expect data to be a data.frame with 2 cols
  #colnames(sim_df) = c('return','t')
  #has 0:t rows ie t+1 rows
  
  
  # init_params is constrainted version (not the ones \in (-inf,inf))
  #NE-STD ones
  init_params[1] = init_params[1] - init_params[3]^-1
  init_params[4] = init_params[4] - init_params[6]^-1
  init_params[1:6] = (init_params[1:6])
  #w
  init_params[7] = (init_params[7])
  #b  (do b first cuz b dependent a constrained)
  init_params[9] =  (-init_params[9]^2+init_params[9]-1/4)/
    (init_params[9]*(init_params[9]-1)*(1-init_params[8]))
  init_params[9] = sqrt(init_params[9])
  #a       a_un = (-a^2+a-1/4)/(a*(a-1))
  init_params[8] =  (-init_params[8]^2+init_params[8]-1/4)/
    (init_params[8]*(init_params[8]-1))
  init_params[8] = sqrt(init_params[8])
  #sigma0
  init_params[10] = (init_params[10])
  
  # result$par = params is in order alpha1,gam1,beta1,alpha2,gam2,beta2,
  #                    pre_w,pre_a,pre_b,sigma0
  print(init_params)
  Results = optim(par = init_params,
                  fn = NE_WGARCH_loglik,#gr = NE_WGARCH_partial_loglik,
                  r0 = data$return[1], r = data$return[-1],
                  control = list(fnscale = -1, maxit = maxit),...
  )
  # NE_WT params
  params = Results$par
  params_1 = abs(params[1:3])
  params_1[1] = params_1[1] + params_1[3]^(-1)
  #inequality \alpha\beta -1 \geq 0
  #  \implies \alpha \geq 1/\beta
  
  params_2 = abs(params[4:6])
  params_2[1] = params_2[1] + params_2[3]^(-1)
  # GARCH params
  w = abs(params[7])
  a = 0.5*(params[8]/sqrt(1+params[8]^2) +1 )
  b = (1-a)*0.5*(params[9]/sqrt(1+params[9]^2) +1 )
  sigma0 = abs(params[10])
  
  Results$par = c(params_1,params_2,w,a,b,sigma0)
  return(Results)
}



