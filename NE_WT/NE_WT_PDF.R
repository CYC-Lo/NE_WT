#PDF of NE-WT

source('NE_W_PDF.R')
NE_WT_PDF = function(x, theta,
                     alpha1,gam1,beta1,
                     alpha2,gam2,beta2){
  #theta is probability negative 
  #param_1 => parameters for >=0
  #param_2 => parameters for <0
  x_idx_neg = which(x < 0)
  x_idx_pos = which(x >=0)
  M = x
  
  if (length(x_idx_neg) != 0) {
    M[x_idx_neg] = theta*NE_W_PDF(-x[x_idx_neg], alpha2,gam2,beta2)
    # print(theta*NE_W_PDF(-x[x_idx_neg], alpha2,gam2,beta2))
  }
  if (length(x_idx_pos) != 0) {
    M[x_idx_pos] = (1-theta)*NE_W_PDF(x[x_idx_pos], alpha1,gam1,beta1)
    # print((1-theta)*NE_W_PDF(x[x_idx_pos], alpha1,gam1,beta1))
  }
  return(M)
}
