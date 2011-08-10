PP_GPCMwle <-
function(u,sl,s,theta0=0,exac=0.001,osch=44,...)
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

#1 -------------------------------------------
P <- function(k,a,th,bet)
{
Zo <- a*(th - bet[1:k])
NeN <- sum(sapply(1:length(bet),function(x) exp(sum(a*(th-bet[1:x])))))

exp(sum(Zo)) / NeN
}#


P1 <- function(k,a,th,bet)
{
# k = answered which category (from 1 to max)
cc <- 1:length(bet)
a * P(k,a,th,bet) * (k - sum(cc * sapply(cc,P,a,th,bet)))
}#  

P2 <- function(k,a,th,bet)
{
cc <- 1:length(bet)
a * P1(k,a,th,bet) * (k - sum(cc * sapply(cc,P,a,th,bet))) - a * P(k,a,th,bet) * (sum(cc * sapply(cc,P1,a,th,bet)))
}#

  
L1 <- function(k,a,th,bet)
        {
        cc <- 1:length(bet)
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
        cc <- 1:length(bet)
        a * sum(cc * sapply(cc,P1,a,th,bet))
        }#


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
cc <- 1:length(bet)
k - sum(cc * sapply(cc,P,a,th,bet))
}#

# 3.Term des Bias
B <- function(k,a,th,bet)
{
cc <- 1:length(bet)
k^2 - 2*k*sum(cc * sapply(cc,P,a,th,bet)) + 2*(sum(cc * sapply(cc,P,a,th,bet)))^2 - sum(cc^2 * sapply(cc,P,a,th,bet))
}#

WT <- function(k,a,th,bet)
{
F(k,a,th,bet) * A(k,a,th,bet) * B(k,a,th,bet)
}# (?)
    
        
GWEIGHT <- function(AA,th,schwellen)
              {
              # AA=slope-vector, th=theta, schwellen = matrix mit schwellenschwierigkeiten je item
              
                gesW <- sapply(1:dim(schwellen)[2],function(schw)				# Geht die Items durch = ?u?ere Summe
              				{
              # vorsicht: hier muss noch das "a" wechseln!
              				kat.W <- sapply(1:length(schwellen[,schw]),function(x) 	# k bis ANZ-Kategorien
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
        cc <- 1:length(bet)
        
        -sum(cc * sapply(cc,P1,a,th,bet))
        }

        
B1 <- function(k,a,th,bet)
        {
        cc <- 1:length(bet)
        
        -2*k * sum(cc * sapply(cc,P1,a,th,bet)) + 4*sum(cc * sapply(cc,P,a,th,bet))*sum(cc * sapply(cc,P1,a,th,bet)) - sum(cc^2 * sapply(cc,P1,a,th,bet))
        }        
        

II1 <- function(k,a,th,bet)
          {
          cc <- 1:length(bet)
          
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
            
            				kat.W <- sapply(1:length(schwellen[,schw]),function(x)
            							{
            							WT_1(k=x,a=AA[schw],th,bet=schwellen[,schw])
            							})
            							
            				Weight <- sum(kat.W)
            				Weight
            	
            			})
            	gesW1
            	}        
        
  
# fisher scoring
count <- 0
ergv <- rep(NA,osch)
      repeat
        {
        delo <- sum(GESL1(slope=sl,kat=u,th=theta0,schwellen=swe)) + sum(GWEIGHT(AA=sl,th=theta0,schwellen=swe)) /(2*sum(GI(AA=sl,th=theta0,schwellen=swe)))
    		delu <- sum(GI(AA=sl,th=theta0,schwellen=swe)) + (sum(GWEIGHT1(AA=sl,th=theta0,schwellen=swe)) / ((2*sum(GI(AA=sl,th=theta0,schwellen=swe)))^2))
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
			

SE <- 1/sqrt(sum(GI(AA=sl,th=theta0,schwellen=swe)))

REP <- list("resp"=u,"estimate"=theta0,"iterations"=count,"estproc"=ergv,"thresh"=swe,"SE"=SE)
class(REP) <- "PPp"
REP
}

