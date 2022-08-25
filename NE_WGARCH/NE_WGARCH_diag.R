#Get diagnostics of the fitted NE_WGARCH

NE_WGARCH_sig = function(MLE_results,data){
  #MLE_results has the mle parameters (constrained form)
  params = MLE_results$par
  alpha1 = params[1]
  gam1 = params[2]
  beta1 = params[3]
  alpha2 = params[4]
  gam2 = params[5]
  beta2 = params[6]
  w = params[7]
  a = params[8]
  b = params[9]
  sigma0 = params[10]

  t = nrow(data)-1
  r = data$return[-1]
  
  #calculate distribution variance
  # sigma_sq_ne = moment_func(2,NE_STD_PDF, c(-Inf, Inf),
  #                           alpha1 = alpha1, gam1 = gam1, beta1 = beta1,
  #                           alpha2 = alpha2, gam2 = gam2, beta2 = beta2)[[1]]
  sigma_sq_ne = 1
  #transform from standard GARCH parameter
  #w^{GARCH} = w^{var \neq 1} \sigma_sq_{ne}
  w =w/sigma_sq_ne
  a =a/sigma_sq_ne
  sigma_sq = c()
  sigma_sq[1] = sigma0^2/sigma_sq_ne
  
  for (i in 2:(t+1) ) {
    sigma_sq[i] = w + a*r[i-1]^2 + b*sigma_sq[i-1]
    if (i==2) {
      print(paste0('w=',w))
      print(paste0('a=',a))
      print(paste0('r[1]=',r[i-1]))
      print(paste0('b=',b))
      print(paste0('sigma_sq[1]=',sigma_sq[i-1]))
    }
  }
  return(sigma_sq)
}

