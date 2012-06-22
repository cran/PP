P <- function(k,a,th,bet)
{
  # bet: threshold parameters
  # k = answered in the kth category
  # this function is applied to each item separately
  bet_s <- as.vector(na.omit(bet))
  Zo    <- a*(th - bet_s[1:k])
  NeN   <- sum(sapply(1:length(bet_s),function(x) exp(sum(a*(th-bet_s[1:x])))))
  
  exp(sum(Zo)) / NeN
}


P1 <- function(k,a,th,bet)
{
  # k = answered which category (from 1 to max)
  bet_s <- as.vector(na.omit(bet))
  cc    <- 1:length(bet_s)
  a * P(k,a,th,bet_s) * (k - sum(cc * sapply(cc,P,a,th,bet_s)))
}  


L1 <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  k*a - a * sum(cc * sapply(cc,P,a,th,bet_s))
}


GESL1  <- function(slope,kat,th,schwellen)
{   # slope = slope VECTOR!
  # kat = xij
  
  gesl1 <- sapply(1:dim(schwellen)[2],function(xx)
  {
    L1(k=kat[xx],a=slope[xx],th,bet=schwellen[,xx])
  })
  gesl1
} 

#-------------------------------------------------------------


II <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  a * sum(cc * sapply(cc,P1,a,th,bet_s))
}


GI <- function(AA,th,schwellen)
{# AA == slope 
  gei <- sapply(1:dim(schwellen)[2],function(schw)
  {
    sona <- as.vector(na.omit(schwellen[,schw])) 
    kat.Inf <- sapply(1:length(sona),function(x) # goes through every threshold
    {
      II(k=x,a=AA[schw],th,bet=schwellen[,schw]) * P(k=x,a=AA[schw],th,bet=schwellen[,schw])
    })
    
    Iteminf <- sum(kat.Inf)
    Iteminf
  })
  gei  	
}
