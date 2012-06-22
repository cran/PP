#1 -------------------------------------------
P <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  Zo <- a*(th - bet_s[1:k])
  NeN <- sum(sapply(1:length(bet_s),function(x) exp(sum(a*(th-bet_s[1:x])))))
  
  exp(sum(Zo)) / NeN
}#


P1 <- function(k,a,th,bet)
{
  # k = choose which category (from 1 to max)
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  a * P(k,a,th,bet) * (k - sum(cc * sapply(cc,P,a,th,bet)))
}#  


P2 <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  a * P1(k,a,th,bet) * (k - sum(cc * sapply(cc,P,a,th,bet))) - a * P(k,a,th,bet) * (sum(cc * sapply(cc,P1,a,th,bet)))
}#


L1 <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  k*a - a * sum(cc * sapply(cc,P,a,th,bet))
}#


GESL1  <- function(slope,kat,th,schwellen)
{   # slope = slope VECTOR!
  # kat = xij
  
  gesl1 <- sapply(1:dim(schwellen)[2],function(xx)
  {
    L1(k=kat[xx],a=slope[xx],th,bet=schwellen[,xx])
  })
  gesl1
} #
#-------------------------------------------------------------

#2 -----------------------------------------------------------
II <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  a * sum(cc * sapply(cc,P1,a,th,bet))
}#


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
}#
#-------------------------------------------------------------

#3 -----------------------------------------------------------  

F <- function(k,a,th,bet)
{
  a^3 * P(k,a,th,bet)
}#

# 2.Term des Bias
A <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc    <- 1:length(bet_s)
  k - sum(cc * sapply(cc,P,a,th,bet))
}#

# 3.Term des Bias
B <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  k^2 - 2*k*sum(cc * sapply(cc,P,a,th,bet)) + 2*(sum(cc * sapply(cc,P,a,th,bet)))^2 - sum(cc^2 * sapply(cc,P,a,th,bet))
}#

WT <- function(k,a,th,bet)
{
  F(k,a,th,bet) * A(k,a,th,bet) * B(k,a,th,bet)
}#


GWEIGHT <- function(AA,th,schwellen)
{
  # AA=slope-vector, th=theta, schwellen = matrix mit schwellenschwierigkeiten je item
  
  gesW <- sapply(1:dim(schwellen)[2],function(schw)				# Geht die Items durch = ?u?ere Summe
  {
    # 
    sona <- as.vector(na.omit(schwellen[,schw])) ## NEU
    kat.W <- sapply(1:length(sona),function(x) 	# k bis ANZ-Kategorien
    {								# innere Summe
      WT(k=x,a=AA[schw],th,bet=schwellen[,schw])
    })
    
    Weight <- sum(kat.W)
    Weight
    
  })
  gesW
}

#-------------------------------------------------------------

#4 -----------------------------------------------------------                      

F1 <- function(k,a,th,bet)
{
  a^3  * P1(k,a,th,bet) 
}


A1 <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc    <- 1:length(bet_s)
  
  -sum(cc * sapply(cc,P1,a,th,bet))
}


B1 <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  
  -2*k * sum(cc * sapply(cc,P1,a,th,bet)) + 4*sum(cc * sapply(cc,P,a,th,bet))*sum(cc * sapply(cc,P1,a,th,bet)) - sum(cc^2 * sapply(cc,P1,a,th,bet))
}        


II1 <- function(k,a,th,bet)
{
  bet_s <- as.vector(na.omit(bet))
  cc <- 1:length(bet_s)
  a * sum(cc * sapply(cc,P2,a,th,bet))
}        


WT_1 <- function(k,a,th,bet)
{
  (A(k,a,th,bet)*B(k,a,th,bet)*2*II(k,a,th,bet)*F1(k,a,th,bet)
   + F(k,a,th,bet)*B(k,a,th,bet)*2*II(k,a,th,bet)*A1(k,a,th,bet)
   + F(k,a,th,bet)*A(k,a,th,bet)*2*II(k,a,th,bet)*B1(k,a,th,bet)
   - F(k,a,th,bet)*A(k,a,th,bet)*B(k,a,th,bet)*2*II1(k,a,th,bet))
}        



GWEIGHT1 <- function(AA,th,schwellen)
{
  
  gesW1 <- sapply(1:dim(schwellen)[2],function(schw)
  {
    sona <- as.vector(na.omit(schwellen[,schw])) ## NEU
    
    kat.W <- sapply(1:length(sona),function(x)
    {
      WT_1(k=x,a=AA[schw],th,bet=schwellen[,schw])
    })
    
    Weight <- sum(kat.W)
    Weight
    
  })
  gesW1
}        
