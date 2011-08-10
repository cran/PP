PP_RMwle <-
function(s,...)
#
{
resm <- ifelse(upper.tri(matrix(1,length(s)+1,length(s)),diag=T),0,1)
a <- rep(1,length(s))
i <- rep(0,length(s))
erg1 <- apply(resm,1,function(o) PP_3PLwle(u=o,a=a,s=s,i=i,...))  

est  <- sapply(erg1,function(x)x$estimate)
SEl  <- sapply(erg1,function(x)x$SE)
ergL <- cbind("Rawscore"=rowSums(resm),"ml_estimate"=est,"SE"=SEl)
backout <- list("result"=erg1,"estimates"=ergL)

class(backout) <- "PPr"
backout
}

