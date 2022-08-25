#CDF of NE-Wei
#Half side of NE-WT


NE_W_CDF = function(x, alpha, gam, beta){
  H_beta = (1-exp(-gam*x^alpha))^beta
  M = 1- (1-H_beta)/exp(H_beta)
  return(M)
}
