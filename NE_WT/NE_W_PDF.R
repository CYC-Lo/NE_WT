#PDF of NE-Wei
#Half side of NE-WT


NE_W_PDF = function(x, alpha, gam, beta){
  exp_gxa = exp(-gam*x^alpha)
  H = (1-exp_gxa)
  #Do exp_gxa first since exp_gxa \approx 0 when x large
  #   this makes x^(alpha-1) huge and overloads the numerics making it Na before
  #   multiplying by 0
  
  #Actually from asymptotic analysis, we know the exp_gxa term goes to zero
  #   while x^(alpha-1) goes to infinity
  #   But exp trumps polynomials \implies if exp_gxa = 0 then m = 0
  #   \implies we make a special case for 0
  m = c()
  if (length(which(exp_gxa == 0))!= 0) {
    m[which(exp_gxa == 0)] = 0
  }
  
  
  if (length(x^alpha == 0) != 0){
    exp_gxa[which(x^alpha ==0 )] = 1
  }
  if (length(which(exp_gxa == 1))!= 0) {
    m[which(exp_gxa == 1)] = 0
    #since exp_gxa = 1, H = 0 and exponential trumps x^(alpha-1)
  }
  
  
  
  idx_norm =  which(exp_gxa != 0 & exp_gxa != 1)
  if (length(idx_norm)!= 0){
  m[idx_norm] = exp_gxa[idx_norm]*(2- H[idx_norm]^beta)*
                alpha*beta*gam*x[idx_norm]^(alpha-1)*
                H[idx_norm]^(beta-1)/exp(H[idx_norm]^beta)
  }
  return(m)  
}
