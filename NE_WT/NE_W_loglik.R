# Log Likelihood equation for NE_W

#Due to optim requirements, it takes all parameters(vector) as first argument




NE_W_logLik = function(params, x){
  alpha = params[1]
  gam = params[2]
  beta = params[3]
  
  H = 1 - exp(-gam*x^alpha)
  # print(paste0( 'params =  ', params) )
  # h = gam*alpha*x^(alpha-1)*exp(-gam*x^alpha)
  
  #Doing log_h to protect numerics from h=0 \implies log(h) = -Inf
  log_h = log(gam)+log(alpha) + (alpha-1)*log(x) - gam*x^alpha
  
  
  
  #Taylor expand log(H) if exp(-gam*x^alpha) too close to 1
  exp_term = exp(-gam*x^alpha)
  log_H = c()
  if(length(which(exp_term !=1)) != 0){
    log_H[which(exp_term != 1)] = log(H[which(exp_term != 1)])
  }
  
  idx_small = which(exp_term == 1) 
  if ( length(idx_small) != 0) {
    log_H[idx_small] = log(gam) + alpha*log(x[idx_small])
  }
  #log(H) = log(1-exp(-gam*x^alpha))
  # since exp(-gam*x^alpha) = 1 - gam*x^alpha + o(gam*x^alpha)
  # \implies  log(H) = log(gam*x^alpha + o(gam*x^alpha))
  #                  \approx log(gam*x^alpha)
  #                 = log(gam) + alpha*log(x)
  
  
  
  # print(paste0( 'H = ', head(H)  ) )
  # print(paste0( 'h = ', head(h)  ) )
  # print(paste0( 'log beta = ', log(beta)  ) )
  # print(paste0( 'length x = ', length(x)  ) )
  # print(paste0( 'which(log_h == -Inf) ',  which(log_h== -Inf)   ))
  # print(paste0( 'which(log_H == -Inf) ',  which(log_H== -Inf)   ))
  # print(paste0('log_h[195] = ', log_h[195]  ))
  # print(paste0('exp_term[1839] = ', exp_term[1839]  ))
  # print(paste0('exp(-gam*x[1839]^(alpha)) = ', exp(-gam*x[1839]^(alpha))  ))
  # print(paste0('-gam*x[195]^(alpha) = ', -gam*x[195]^(alpha)  ))
  # print('log_h len')
  # print(length(log_h))
  # print('log_H len')
  # print(length(log_H))
  L = length(x)*log(beta)+sum( log_h + (beta -1)*log_H + log(2-H^beta) - H^beta)
  return(L)
}
#Bug spot

NE_W_alpha_partial_loglik = function(params,x){
  alpha = params[1]
  gam = params[2]
  beta = params[3]
  
  H = 1- exp(-gam*x^alpha)
  sum_1 = sum(
    log(x)*(1+ (H^-1 -2)* gam *x^alpha ) 
    )
  sum_2 =sum(
    (H^-1 -1)*gam*x^alpha *log(x)*( (2-H^beta) - 2/(2-H^beta) )  
    )
  L = length(x)/alpha + sum_1 + beta*sum_2
  return(L)
}

NE_W_gam_partial_loglik = function(params,x){
  alpha = params[1]
  gam = params[2]
  beta = params[3]
  
  H = 1- exp(-gam*x^alpha)
  sum_1 = sum(
    (H^-1 -2)*x^alpha
    )
  sum_2 = sum(
    (H^-1 -1)*x^alpha * ( (2-H^beta) - 2/(2-H^beta) )  
    )
  L = length(x)/gam + sum_1 + beta*sum_2
  return(L)
}


NE_W_beta_partial_loglik = function(params,x){
  alpha = params[1]
  gam = params[2]
  beta = params[3]
  
  H = 1- exp(-gam*x^alpha)
  L = length(x)/beta +
    sum(
      log(H)*( (2-H^beta) - 2/(2-H^beta) )
      )
  return(L)
}

NE_W_partial_loglik = function(params,x){
  return(
    c(NE_W_alpha_partial_loglik(params,x),
      NE_W_gam_partial_loglik(params,x),
      NE_W_beta_partial_loglik(params,x)
      )
    )
  
}


NE_W_MLE = function(data, init_params = c(2,2,2),...){
  return(Results = optim(par = init_params, fn = NE_W_logLik,
                         gr = NE_W_partial_loglik, x = data,
                         control = list(fnscale = -1)
                         ) 
         )
}
