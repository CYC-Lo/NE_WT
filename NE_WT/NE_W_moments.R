#Moment integrator for NE_Weibull

source('NE_W_PDF.R')
moment_func = function(r, pdf, domain = c(-Inf,Inf),...){
  # r = rth moment
  # pdf is the pdf of the distribution of interest
  # domain is the domain for the PDF
  new_func = function(x,...){
    mu = x^r * pdf(x,...)
    return(mu)
  }
  moment_val = integrate(new_func, lower = domain[1], upper = domain[2],...)
  return(moment_val)
}


NE_W_moment = function(r,domain = c(0,Inf),...){
  moment = moment_func(r,pdf = NE_W_PDF, domain = domain, ...)
  return(moment)
}
  