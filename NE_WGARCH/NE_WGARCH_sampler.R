# NE_WGARCH sampler

setwd('~/Documents/MSc_thesis/NE_WGARCH/')

source('../NE_WT/NE_STD_sampler.R',chdir = TRUE)
source('../NE_WT/NE_STD_PDF.R',chdir = TRUE)


NE_WGARCH_sampler = function(t,
                          alpha1,gam1,beta1,
                          alpha2,gam2,beta2,
                          w,a,b,sigma0,
                          error = 1e-6,init_error = 0.05){
  # t is number of time points
  #NE_STD param = alpha1,gam1,beta1,alpha2,gam2,beta2 
  
  #calculate distribution variance
  sigma_sq_ne = moment_func(2,NE_STD_PDF, c(-Inf, Inf),
                         alpha1 = alpha1, gam1 = gam1, beta1 = beta1,
                         alpha2 = alpha2, gam2 = gam2, beta2 = beta2)[[1]]
  #sample from NE_STD with variance sigma_ne
  x = NE_STD_sampler(t+1, alpha1,gam1,beta1,
                     alpha2,gam2,beta2)[,1]
  
  #transform from standard GARCH parameter
  #w^{GARCH} = w^{var \neq 1} \sigma_sq_{ne}
  w =w/sigma_sq_ne
  a =a/sigma_sq_ne
  
  #set up r and sigma, they go  from 1:(t+1). time index is 0:t
  #transform into sigma_0^2
  sigma_sq = c()
  sigma_sq[1] = sigma0^2/sigma_sq_ne
  r = c()
  r[1] = sqrt(sigma_sq[1])*x[1]
  # now we have r0, sigma_sq0
  
  for (i in 2:(t+1) ) {
    sigma_sq[i] = w + a*r[i-1]^2 + b*sigma_sq[i-1]
    r[i] = sqrt(sigma_sq[i])*x[i]

  }
  interm = list(w_garch = w, a_garch = a, sigma_sq_ne = sigma_sq_ne)
  input = as.list(match.call())
  return(list(return = r, var = sigma_sq, t = 0:t, x = x,input = input,interm = interm))
  
}
