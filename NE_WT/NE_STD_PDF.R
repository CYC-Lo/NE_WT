#PDF of NE-STD

source('NE_WT_PDF.R')
NE_STD_PDF = function(x,
                      alpha1,gam1,beta1,
                      alpha2,gam2,beta2){
  #param_1 => parameters for >=0
  #param_2 => parameters for <0
  
  
  #theta is probability negative, decided by E[X_1]/(E[X_1]- E[X_2]) 
  #Calculate theta
  mean_1 = NE_W_moment(1,c(0,Inf),
                       alpha = alpha1, gam = gam1,
                       beta = beta1)[[1]]
  mean_2 = -NE_W_moment(1,c(0,Inf),
                        alpha = alpha2, gam = gam2,
                        beta = beta2)[[1]]
  theta = mean_1/(mean_1-mean_2)
  
  return(NE_WT_PDF(x, theta,
                   alpha1,gam1,beta1,
                   alpha2,gam2,beta2)
  )
  
}
