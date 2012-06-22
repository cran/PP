PP_3PLmle <-
function(u,a,s,i,theta0=0,exac=0.001,osch=22, ...)
{
#u - response vector - numeric - binary format is necessary (0,1) 
#a - vector of slope parameters - associated with the items used
#s - vector of difficulty parameters
#i - vector with lower asymptote parameters
#theta0 - starting number for parameter estimation
#exac - accuracy - at which difference between two consecutive estimations the algorithm should stop

#ctrl
co_3PL(u,a,s,i,theta0)


count <- 0 
ergv <- rep(NA,osch)
if(all(u==1))
  {
  theta0 <- Inf  
  SE <- NA
  } else if(all(u==0))
      {
      theta0 <- -Inf 
      SE <- NA
      } else 
          {
          repeat
            {
            del <- DELTA_mle(a,s,i,th=theta0,u)$deLt
            if(is.na(del))
              {
  						theta0 <- NA
              SE     <- NA  
  						break
              }
            
            if(abs(del)>2){del <- del/abs(del) * 2}
            thetan <- theta0 + del
            theta0 <- thetan
            count <- count+1
            ergv[count] <- theta0
            if(abs(del) < exac | count >= osch)
              {
              SE <- 1/sqrt(DELTA_mle(a,s,i,th=theta0,u)$InF)
              break
              }
            } 
          }


REP <- list("resp"=u,"estimate"=theta0,"iterations"=count,"estproc"=ergv,"SE"=SE)
class(REP) <- c("PPd")

REP
}

