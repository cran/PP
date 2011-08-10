summary.PPd <-
function(object,...)
{
  
cat("Response vector: \n")  
print(object$resp)
cat("\nNumber of Items solved: \n")  
print(sum(object$resp))
cat("\nEstimate of person ability and SE: \n")
print(data.frame("PP"=object$estimate,"SE"=object$SE,row.names=""))
cat("\nNumber of iterations: \n")
print(object$iterations)
cat("\nEstimation progress: \n")

prog <- object$estproc[-which(is.na(object$estproc))]
print(prog)
  
}

