\name{PP_3PLmle}
\alias{PP_3PLmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of ML person parameters for the 3-PL model
}
\description{
This function estimates person parameters for a given response vector and fixed item parameters (slope, difficulty, lower-asymptote) of the 3-PL model under a maximum-likelihood approach.
}
\usage{
PP_3PLmle(u, a, s, i, theta0 = 0, exac = 0.001, osch = 22, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
numeric response vector - binary format is required (0,1)
}
  \item{a}{
numeric vector of slope parameters
}
  \item{s}{
numeric vector of difficulty parameters (item difficulty)
}
  \item{i}{
numeric vector with lower-asymptote parameters
}
  \item{theta0}{
  starting value for the person parameter (default = 0)
}
  \item{exac}{
accuracy criterion - at which difference between two consecutive estimations the algorithm should stop
}
  \item{osch}{
maximum number of iterations (default is arbitrarily set to 22)
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}


\value{
Returns an object of class "PPd" which contains the following values:

\item{resp}{entered response vector}
\item{estimate}{person parameter estimate}
\item{iterations}{number of iterations}
\item{estproc}{estimate for every iteration step - the resulting vector is of length(1:osch) - usually the algorithm stops before reaching the iteration limit, so there are NA's resulting.}
\item{SE}{standard error}

}
\references{
Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.
}
\author{
Manuel Reif
}


\seealso{
\code{\link{PP_3PLwle}}
}
\examples{
PP_3PLmle(u=c(1,0,1,0),a=c(1,1,1,1),s=c(-2,-1.23,1.11,3.48),i=rep(0,4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
