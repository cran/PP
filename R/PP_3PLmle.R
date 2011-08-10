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
#1
lelist <- sapply(list(u,a,s,i),length)
if(any(lelist != lelist[1])){stop("Check length of input vectors \n")}
#2 
if(length(theta0) != 1){stop("theta0 must be of length 1 \n")}
#3
if(any(!(u %in% c(0,1)))){stop("u is not of format (0,1) as required \n")}

  
DELTA_mle <- function(a,s,i,th,u) 
    {
		# 3PL model 
		Zae <- expression(exp(a*(th-s)))
		P <- i+(1-i) * eval(Zae) / (1+eval(Zae))
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

