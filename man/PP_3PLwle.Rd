\name{PP_3PLwle}
\alias{PP_3PLwle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of WL Person Parameters for the 3-PL model
}
\description{
This function estimates person parameters for a given response vector and fixed parameters (slope, difficulty, lower-asymptote) of the 3-PL model under Warms weighted maximum-likelihood approach (Warm 1989).
}
\usage{
PP_3PLwle(u, a, s, i, theta0 = 0, exac = 0.001, osch = 33, ...)
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
\details{
%%  ~~ If necessary, more details than the description above ~~
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

Warm, Thomas A. (1989). Weighted Likelihood Estimation Of Ability In Item Response Theory. Psychometrika, 54, 427-450.
}
\author{
Manuel Reif
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PP_3PLwle(u=c(1,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
