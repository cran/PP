co_3PL <- 
function(u,a,s,i,theta0)
  {
    #1
    lelist <- sapply(list(u,a,s,i),length)
    if(any(lelist != lelist[1])){stop("Check length of input vectors \n")}
    
    #2 
    if(length(theta0) != 1){stop("theta0 must be of length 1 \n")}
    
    #3
    if(any(!(u %in% c(0,1)))){stop("u is not of format (0,1) as required \n")}
      
  }