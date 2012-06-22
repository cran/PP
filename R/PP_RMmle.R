PP_RMmle <-
function(s,expol=FALSE,...)
#
{
#
if(expol & length(s) <= 2){stop("Cannot perform extrapolation with #items <= 2")}
  
resm <- ifelse(upper.tri(matrix(1,length(s)+1,length(s)),diag=T),0,1)
a <- rep(1,length(s))
i <- rep(0,length(s))
erg1 <- apply(resm,1,function(o) PP_3PLmle(u=o,a=a,s=s,i=i,...))

wexpol <- "no expol performed"

# EXPOL
if(expol)
  {
  est  <- sapply(erg1,function(x)x$estimate)
  ergL <- cbind("Rawscore"=rowSums(resm),"ml_estimate"=est)
  ergLk <- ergL[-c(1,nrow(ergL)),]
  
  if((length(est)-2) >= 4)
    {
    rmi <- interpSpline(ergLk[,1],ergLk[,2])
    pre <- predict(rmi,c(0,max(ergLk[,1])+1))
    ergL[c(1,nrow(ergL)),2] <- pre$y
    wexpol <- "performed spline expol "
    SEl  <- sapply(erg1,function(x)x$SE)
    ergL <- cbind(ergL,"SE"=SEl)
    } else
        {
        rmi <- lm(ml_estimate ~ Rawscore,data=data.frame(ergLk))
        nd <- data.frame(Rawscore=c(0,(max(ergLk[,1])+1)))
        pre <- predict(rmi,newdata=nd)
        ergL[c(1,nrow(ergL)),2] <- pre
        wexpol <- "performed linmod expol"
        SEl  <- sapply(erg1,function(x)x$SE)
        ergL <- cbind(ergL,"SE"=SEl)
        }
  } else
      {
      est  <- sapply(erg1,function(x)x$estimate)
      SEl  <- sapply(erg1,function(x)x$SE)
      ergL <- cbind("Rawscore"=rowSums(resm),"ml_estimate"=est,"SE"=SEl)  
      }


backout <- list("result"=erg1,"estimates"=ergL,"extrap"=wexpol)

class(backout) <- "PPr"
backout
}

