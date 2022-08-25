#Moment integrator for NE_WT

source('NE_W_moments.R')
NE_WT_moment = function(r,domain = c(-Inf,Inf),...){
  moment = moment_func(r, pdf = NE_WT_PDF, domain = domain,...)
  return(moment)
}

