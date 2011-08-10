PP_3PLwle <-
function(u,a,s,i,theta0=0,exac=0.001,osch=33, ...)
{
  
#ctrl
#1
lelist <- sapply(list(u,a,s,i),length)
if(any(lelist != lelist[1])){stop("Check length of input vectors \n")}
#2 
if(length(theta0) != 1){stop("theta0 must be of length 1 \n")}
#3
if(any(!(u %in% c(0,1)))){stop("u is not of format (0,1) as required \n")}
  
  
DELTA_wle <- function(a,s,i,th,u) 
    {
		# th = theta (ability parameter - in the first run a starting value will be assigned to th)
		# u is the binary response vector
		
		# 3PL model 
		Zae <- expression(exp(a*(th-s)))
		P <- i+(1-i) * eval(Zae) / (1+eval(Zae))
		# 1st derivate of P
		P.1ab <- a*(1-P)*(P-i)/(1-i)
		# 2nd derivate of P
		P.2ab <- (a/(1-i)) * (P.1ab - 2*P.1ab*P + P.1ab*i)
		# 3rd derivate of P
		P.3ab <- - a*((i-2*P +1) * P.2ab - 2*P.1ab^2) / (i-1)
		# shortcuts
		Emp <- 1-P
		Pmi <- P-i
		Emi <- 1-i
		# 1st derivate of the logL
		l1 <- sum(a*P*Emp*(u-P)/(P*Emp)*(Pmi/Emi/P)) 
		# I
		 I <- sum(P.1ab^2 / (P*Emp))
		# J
		 J <- sum((P.1ab*P.2ab) / (P*Emp))
		#
		 Te1.I1 <- 2*P.1ab*P.2ab/(P*Emp) 
		 Te2.I1 <- (P.1ab^2 * (P.1ab*Emp - P*P.1ab)) / (P*Emp)^2
		 I.1ab <- sum(Te1.I1 - Te2.I1)
		#
		J.1ab <- sum(-((P-1)*P*P.2ab^2 - (P-1)*P*P.3ab*P.1ab + (2*P-1)*P.1ab^2*P.2ab)/(P^2 * (P-1)^2))
		# put it all together
		O <- (l1 + (J/(2*I)))
		U <- I + (I*J.1ab - I.1ab*J) / (2*I^2)
		delta <- O/U
		list(deLt=delta,InF=I)
		}  
  

count <- 0
ergv <- rep(NA,osch)
repeat
	{
	del <- DELTA_wle(a,s,i,th=theta0,u)$deLt
  if(is.na(del))
    {
  	theta0 <- NA
    SE     <- NA  
		break
    }
	if(abs(del)>2){del <- del/abs(del) * 2}
  thetan <- theta0 + del
	theta0 <- thetan
  count <- count + 1
  ergv[count] <- theta0
	if(abs(del) <= exac | count >= osch)
    {
    SE <- 1/sqrt(DELTA_wle(a,s,i,th=theta0,u)$InF)
    break
    }
	
	}

REP <- list("resp"=u,"estimate"=theta0,"iterations"=count,"estproc"=ergv,"SE"=SE)
class(REP) <- c("PPd")
REP
}

