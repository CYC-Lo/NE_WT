#Log Likelihood function for NE-WGARCH
#No partial derivatives at the moment

source('../NE_WT/NE_STD_loglik.R',chdir = TRUE)


NE_WGARCH_loglik = function(params,r,r0){
  # params is in order alpha1,gam1,beta1,alpha2,gam2,beta2,
  #                    pre_w,pre_a,pre_b,sigma0
  
  # Basically NE_WT params then GARCH params
  
  # r is data, r or epsilon from 1:T (not the x from notation)
  # r0 is epsilon 0, the data point before first day of interest
  
  # NE_WT params 
  #Need to transform to ensure all > 0
  #   as well as alpha*beta -1 \geq 0
  params_1 = exp(params[1:3])
  params_1[1] = params_1[1] + params_1[3]^(-1)
  params_2 = exp(params[4:6])
  params_2[1] = params_2[1] + params_2[3]^(-1)
  
  # GARCH params
  #Need to transform to fit constraint
  w = exp(params[7])
  a = 0.5*(params[8]/sqrt(1+params[8]^2) +1 )
  b = (1-a)*0.5*(params[9]/sqrt(1+params[9]^2) +1 )
  sigma0 = exp(params[10])
  
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


NE_WGARCH_partial_loglik = function(params,r,r0){
  # params is in order alpha1,gam1,beta1,alpha2,gam2,beta2,
  #                    pre_w,pre_a,pre_b,sigma0
  
  # NE_WT params
  params_1 = exp(params[1:3])
  params_1[1] = params_1[1] + params_1[3]^(-1)
  params_2 = exp(params[4:6])
  params_2[1] = params_2[1] + params_2[3]^(-1)
  # GARCH params
  w = exp(params[7])
  a = 0.5*(params[8]/sqrt(1+params[8]^2) +1 )
  b = (1-a)*0.5*(params[9]/sqrt(1+params[9]^2) +1 )
  sigma0 = exp(params[10])
  
  T_len = length(r)
  # =======================
  #transform r to x
  # step 1
  sigma_sq = c()
  sigma_sq[1] = w + a*r0^2 + b*sigma0^2
  for (i in 2:T_len) {
    sigma_sq[i] = w + a*r[i-1]^2 + b*sigma_sq[i-1]
  }
  # step 2
  x = r/sqrt(sigma_sq)
  
  #========================
  #split into N_1 and N_2 for L_{N_1} and L_{N_2}
  x_1 = x[x >= 0]
  x_2 = -(x[x < 0])
  # x_1 is postive real line data
  # x_2 is  absolute value of negatice real line data
  N_1 = length(x_1)
  N_2 = length(x_2)
  
  #========================
  #Calculate partial wrt to NE_WT params
  partial_1 = NE_W_partial_loglik(params_1,x_1)
  partial_2= NE_W_partial_loglik(params_2,x_2)
  
  
  #========================
  #calculate dx_dA
  # dx_dA = matrix of column labels(dx_t_dw,dx_t_da,dx_t_db)
  # each row j\in[1,T] is dx_j
  # treated as Tx3 by r
  
  dx_dA = (-0.5*x/sigma_sq) *
    matrix(c(rep(1,T_len), 
             c(r0^2, r[1:(T_len-1)]^2),
             c(sigma0^2, sigma_sq[1:(T_len-1)])
             ),
           ncol = 3)

  #========================
  #calculate partial wrt to GARCH params
  partial_garch = NE_W_A_partial_loglik(params_1,x_1,dx_dA[x>=0,])+
    NE_W_A_partial_loglik(params_2,x_2,dx_dA[x<0,])-
    0.5/sigma_sq%*% matrix(c(rep(1,T_len),
                             c(r0^2, r[1:(T_len-1)]^2),
                             c(sigma0^2, sigma_sq[1:(T_len-1)])
                             ),
                           ncol = 3)
  
  #========================
  #calculate partial wrt to sigma0
  # decide which sign x1 is
  if(x[1]>=0){
    partial_sigma0 = -b*sigma0/sigma_sq[1]*
      (1+ x[1]*NE_W_A_partial_loglik(params_1,x[1],1))
  } else {
    partial_sigma0 = -b*sigma0/sigma_sq[1]*
      (1-x[1]*NE_W_A_partial_loglik(params_2,-x[1],1))
  }
  
  
  # params is in order alpha1,gam1,beta1,alpha2,gam2,beta2,
  #                    w,a,b,sigma0
  return(c(partial_1,partial_2,
           partial_garch,partial_sigma0)
         )
  
}



NE_W_A_partial_loglik = function(params,x,dx_dA){
  # params are alpha,gam,beta
  alpha = params[1]
  gam = params[2]
  beta = params[3]
  #dx_dA = matrix of column labels(dx_t_dw,dx_t_da,dx_t_db)
  #each row j\in[1,T] is dx_j
  # treated as Tx3 by r
  
  h_H = gam*alpha*x^(alpha-1)/(exp(gam*x^alpha)-1) 
  
  
  sum_1 = (alpha-1 + gam*alpha*x^alpha)/x + h_H
  sum_1 = sum_1 %*% dx_dA
  H_beta = (1 - exp(-gam*x^alpha))^beta 
  sum_2 = h_H*((2-H_beta) - 2/(2-H_beta))
  sum_2 = sum_2 %*% dx_dA
  
  # returns vector of (dL_dw,dL_da,dL_db)
  return(sum_1 + beta*sum_2)
}




NE_WGARCH_MLE = function(data, init_params = c(1.5,1.5,1.5,1.5,1.5,1.5,
                                               0.5,0.4,0.4,0.5), maxit = 1000,...){
  #expect data to be a data.frame with 2 cols
  #colnames(sim_df) = c('return','t')
  #has 0:t rows ie t+1 rows
  
  
  # init_params is constrainted version (not the ones \in (-inf,inf))
  #NE-STD ones
  init_params[1] = init_params[1] - init_params[3]^-1
  init_params[4] = init_params[4] - init_params[6]^-1
  init_params[1:6] = log(init_params[1:6])
  #w
  init_params[7] = log(init_params[7])
  #b  (do b first cuz b dependent a constrained)
  init_params[9] =  (-init_params[9]^2+init_params[9]-1/4)/
                    (init_params[9]*(init_params[9]-1)*(1-init_params[8]))
  init_params[9] = sqrt(init_params[9])
  #a       a_un = (-a^2+a-1/4)/(a*(a-1))
  init_params[8] =  (-init_params[8]^2+init_params[8]-1/4)/
                    (init_params[8]*(init_params[8]-1))
  init_params[8] = sqrt(init_params[8])
  #sigma0
  init_params[10] = log(init_params[10])
  
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
  params_1 = exp(params[1:3])
  params_1[1] = params_1[1] + params_1[3]^(-1)
  #inequality \alpha\beta -1 \geq 0
  #  \implies \alpha \geq 1/\beta
  
  params_2 = exp(params[4:6])
  params_2[1] = params_2[1] + params_2[3]^(-1)
  # GARCH params
  w = exp(params[7])
  a = 0.5*(params[8]/sqrt(1+params[8]^2) +1 )
  b = (1-a)*0.5*(params[9]/sqrt(1+params[9]^2) +1 )
  sigma0 = exp(params[10])

  Results$par = c(params_1,params_2,w,a,b,sigma0)
  return(Results)
}



