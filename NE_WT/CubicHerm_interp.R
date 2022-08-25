



#a node is a vector of (p,u,s) = (x-axis, F(p), f(p))
construct_node = function(p,CDF,PDF,...){
  #p is the actual sample value we want
  #CDF is CDF function
  #PDF is PDF function
  #... takes distributions params
  
  #returns a node
  return(
    t(matrix(c(p,CDF(p,...),PDF(p,...)), ncol = 3))
    )
}
node_lst2mat = function(node_lst){
  # takes node list to node mat
  #node list: each element is a node
  #node mat: each col is a node = (p,F(p),f(p))
  return(matrix(unlist(node_lst),nrow = 3))
}
  
  
error_split = function(n_l,n_r,CDF,PDF,error = 1e-6,init_error = 0.05,...){
  #n_l is the node on left bound
  #n_r is the node on right bound
  #CDF is CDF function
  #PDF is PDF function
  #error is acceptable error(between real and interp)
  #init error is F(p_{i+1}) - F(p_{i})

  #returns list of nodes in between n_l,n_r;
  # excludes n_l & n_r
  if ( (n_r[2]-n_l[2])> init_error) {
    p = (n_r[1]+n_l[1])/2
    n_s = construct_node(p,CDF,PDF,...)
    return(
      c( c(error_split(n_l,n_s,CDF,PDF,error,init_error,...)),
         list(n_s),
         c(error_split(n_s,n_r,CDF,PDF,error,init_error,...))
         )
      )
  } else {
    return( list())
  }
}


mono_split = function(n_l,n_r, CDF , PDF ,error = 1e-6,...){
  #n_l is node on left bound
  #n_r is node on right bound
  #CDF is CDF function
  #PDF is PDF function
  #error is acceptable error(between real and interp)
  
  #returns list of nodes in between n_l,n_r;
  # excludes n_l & n_r
  
  criterion = n_r - n_l
  criterion = 3*criterion[1]/criterion[2]
  p = (n_r[1]+n_l[1])/2
  n_s = construct_node(p, CDF, PDF,...)
  if (n_r[1] - n_l[1] <= error ){
    # termination check
    return(list())
  } else if (n_r[3]^-1 > criterion  | n_l[3]^-1>criterion) {
    # monotonicity check
    # if condi = true => not monotonic

    return(
      c( c(mono_split(n_l,n_s,CDF,PDF,error,...)), 
         list(n_s),
         c(mono_split(n_s,n_r,CDF,PDF,error,...))
      )   
    )
  } else {
    # #begin
    spl = splinefunH(x = c(n_r[1],n_l[1]),
                     y = c(n_r[2],n_l[2]),
                     m = c(n_r[3],n_l[3]))
    if ( abs(spl(n_s[1]) - CDF(n_s[1],...)) > error   ) {
      # acceptable error check
      return(
        c( c(mono_split(n_l,n_s,CDF,PDF,error,...)),
           list(n_s),
           c(mono_split(n_s,n_r,CDF,PDF,error,...))
        )
      )
    } else {
      return( list())
    }
    # #end
    return(list())
  }
}

loop_split = function(node_lst,split_func,CDF,PDF,error = 1e-6,...){
  # loops on node_lst to split each interval according to split_func
  
  final_lst = list(node_lst[[1]])
  for (i in 1: (length(node_lst)-1) ){
    n_l = node_lst[[i]]
    n_r = node_lst[[i+1]]
    # print(paste0('n_l',n_l))
    # print(paste0('n_r',n_r))
    final_lst = c(final_lst,
                  c(split_func(n_l,n_r,CDF,PDF,error,...)),
                  list(n_r)
                  )
    
  }
  return(final_lst)
}

