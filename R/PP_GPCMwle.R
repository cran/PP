PP_GPCMwle <-
  function(u,sl,s,theta0=0,exac=0.001,osch=44,...)
  {
    
    co_gpcm(u=u,sl=sl,s=s,theta0=theta0)
    
    # rebuild s
    swe <- rbind(0,s)  
    colnames(swe) <- paste("I",1:ncol(swe),sep="")
    rownames(swe) <- paste("S",1:nrow(swe),sep="")
    # ---- 

    
    
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
