DELTA_mle <- function(a,s,i,th,u) 
{
  # 3PL model 
  Zae <- exp(a*(th-s))
  
  P <- i+(1-i) *  Zae / (1+ Zae)
  # shortcuts
  Emp <- 1-P
  Pmi <- P-i
  Emi <- 1-i
  Wij <- P*(1-P)
  # put it all together
  O <-  sum(a*Wij*((u-P)/Wij)*(Pmi/Emi/P))
  U <-  sum(a^2*Wij*(Pmi/Emi/P)^2) 
  list(deLt = O/U,InF = U)
}









