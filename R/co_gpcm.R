co_gpcm <- 
  function(u,sl,s,theta0)
  {
   
  #  1
  if(min(u) <= 0){stop(paste("Lowest category must be 1 and not", min(u)))}
  
  #  2
  lelist <- sapply(list(u,sl,s[1,]),length)
  
  if(any(lelist != lelist[1])){stop("Check length/dim of input vectors/matrix \n")}
  
  #  3 
  if(length(theta0) != 1){stop("theta0 must be of length 1 \n")}
  
  #  4
  nulist <- sapply(list(u,sl,s,theta0),function(NUM) !is.numeric(NUM))
  if(any(nulist)){stop("Numeric input ... please!")}
  
  }


