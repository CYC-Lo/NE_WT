#Sampler for NE-Wei

source('NE_W_CDF.R')
source('NE_W_PDF.R')
source('CubicHerm_interp.R')


NE_W_sampler = function(n,
                         alpha,gam,beta,
                         error = 1e-6,init_error = 0.05){
  #n is number of samples wanted
  #error and init error are numerical scheme parameters
  
  #sampler on CDF
  u = runif(n)
  u_sort = sort(u)
  
  n_l = construct_node(0,NE_W_CDF,NE_W_PDF,alpha= alpha,gam= gam,beta= beta)
  
  #Ensures we interpolate instead of extrapolate
  i = 20
  n_r = construct_node(i,NE_W_CDF,NE_W_PDF,alpha= alpha,gam = gam,beta= beta)
  while(n_r[2]<max(u)){
    i = i+1
    n_r = construct_node(i,NE_W_CDF,NE_W_PDF,alpha= alpha,gam= gam,beta= beta)
  }
  
  
  vals = loop_split(c(list(n_l),list(n_r)),error_split,NE_W_CDF,NE_W_PDF,
                    alpha = alpha,gam = gam,beta = beta)
  
  vals = loop_split(vals, mono_split,NE_W_CDF,NE_W_PDF,alpha = alpha, gam = gam,beta = beta)
  vals = node_lst2mat(vals)
  # return(vals)
  
  #To catch out the values too similar which the spline function sees are same
  cut_off_l = rev(which(vals[2,]<=0 + error))[1]
  cut_off_r = which(vals[2,]>=1 - error)[1]
  vals_spline = vals[,cut_off_l:cut_off_r]
  spl = splinefunH(vals_spline[2,],vals_spline[1,],vals_spline[3,]^-1)
  
  p = u
  p_idx_l = which(p <= vals[2,cut_off_l])
  p_idx_r = which(p >= vals[2,cut_off_r])
  #Find spline idx
  if(length(p_idx_l) != 0 & length(p_idx_r) != 0 ){
    p_idx_spline = (1:n)[-c(p_idx_l, p_idx_r)]
  } else if (length(p_idx_l) != 0) {
    p_idx_spline = (1:n)[-c(p_idx_l)]
  } else if (length(p_idx_r) != 0){
    p_idx_spline = (1:n)[-c(p_idx_r)]   
  } else {
    p_idx_spline = (1:n)  
  }
  p[p_idx_spline] = spl(p[p_idx_spline])
  
  #Use linear interpolation suggested by paper where not applicable
  slope_l = (vals[1,cut_off_l] - vals[1,1])/(vals[2,cut_off_l]- vals[2,1])
  p[p_idx_l] = slope_l*p[p_idx_l]
  slope_r = (vals[1,dim(vals)[2]] - vals[1,cut_off_r])/(vals[2,dim(vals)[2]] - vals[2,cut_off_l])
  p[p_idx_r] = slope_r*p[p_idx_r] + 1
  
  return(matrix( c(p,u), ncol = 2 ) )
  # return(spl(u_sort)[match(u_sort,u)]  )
  
}
