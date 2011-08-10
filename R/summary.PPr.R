summary.PPr <-
function(object,...)
{
cat("\n Person Parameter Estimates \n")  
  
if(any(names(object)=="extrap"))
{ 
cat("\n",object$extrap,"\n")  
}
  
cat("\n Number of iterations: \n")  
ite  <- sapply(object$result,function(x)x$iteration)
print(ite)
  

cat("\n Estimates for every possible rawscore: \n")
  
object$estimates
  

}

