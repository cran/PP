PP_GPCMmle <-
function(u,sl,s,theta0=0,exac=0.001,osch=33,...)
{

#ctrl
#1
lelist <- sapply(list(u,sl,s[1,]),length)
if(any(lelist != lelist[1])){stop("Check length/dim of input vectors/matrix \n")}
#2 
if(length(theta0) != 1){stop("theta0 must be of length 1 \n")}  

  
# rebuild s
swe <- rbind(0,s)  
colnames(swe) <- paste("I",1:ncol(swe),sep="")
rownames(swe) <- paste("S",1:nrow(swe),sep="")
# ---- 

# make fuctions to estimate theta
#
P <- function(k,a,th,bet)
{
# bet: threshold parameters
Zo <- a*(th - bet[1:k])
NeN <- sum(sapply(1:length(bet),function(x) exp(sum(a*(th-bet[1:x])))))

exp(sum(Zo)) / NeN
}


P1 <- function(k,a,th,bet)
{
# k = answered which category (from 1 to max)
cc <- 1:length(bet)
a * P(k,a,th,bet) * (k - sum(cc * sapply(cc,P,a,th,bet)))
}  
  
  
L1 <- function(k,a,th,bet)
        {
        cc <- 1:length(bet)
        k*a - a * sum(cc * sapply(cc,P,a,th,bet))
        }
  
  
GESL1  <- function(slope,kat,th,schwellen)
              { 	# slope = slope VECTOR!
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
        cc <- 1:length(bet)
        a * sum(cc * sapply(cc,P1,a,th,bet))
        }


GI <- function(AA,th,schwellen)
        {# AA == slope 
        gei <- sapply(1:dim(schwellen)[2],function(schw) # greps every item
          	{
        
        		kat.Inf <- sapply(1:length(schwellen[,schw]),function(x) # goes through every threshold
        					{
        					II(k=x,a=AA[schw],th,bet=schwellen[,schw]) * P(k=x,a=AA[schw],th,bet=schwellen[,schw])
        					})
        					
        		Iteminf <- sum(kat.Inf)
        		Iteminf
        		})
        gei		
        }
  
  
# fisher scoring
count <- 0
ergv <- rep(NA,osch)
  if(all(u == 1))
		{
		theta0 <- -Inf
		} else if(all(u == dim(swe)[1]))
			{
			theta0 <- Inf
			} else 
			{
      repeat
        {
        delo <- sum(GESL1(slope=sl,kat=u,th=theta0,schwellen=swe))
				delu <- sum(GI(AA=sl,th=theta0,schwellen=swe))
        DEL <- delo/delu
        if(is.na(DEL))
          {
  				theta0 <- NA
          SE     <- NA  
					break
          }
        if(abs(DEL)>2){DEL <- DEL/abs(DEL) * 2}
        
        thNEW  <- theta0 + DEL
        ergv[count] <- thNEW
        theta0 <- thNEW 
        if(abs(DEL) <= exac | count >= osch){break}
        count <- count + 1
        }
			}

SE <- 1/sqrt(sum(GI(AA=sl,th=theta0,schwellen=swe)))

REP <- list("resp"=u,"estimate"=theta0,"iterations"=count,"estproc"=ergv,"thresh"=swe,"SE"=SE)
class(REP) <- "PPp"
REP
}

