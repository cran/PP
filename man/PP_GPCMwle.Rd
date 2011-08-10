\name{PP_GPCMwle}
\alias{PP_GPCMwle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of WL Person Parameters for the GPCM
}
\description{
This function estimates Person Parameters for a given response vector (polychotom) and fixed parameters ("slope" and threshold-parameters) of the GPCM (Muraki 1992) under Warms weighted maximum-likelihood approach.
}
\usage{
PP_GPCMwle(u, sl, s, theta0 = 0, exac = 0.001, osch = 44, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
numeric response vector - contains the responses of the examinee (typically integer values; an answer to the lowest category must be denoted as 1)
}
  \item{sl}{
numeric vector of slope parameters (one per item - these are assumed to be constant across all thresholds of one item)
}
  \item{s}{
matrix of threshold parameters (e.g.: for an item with 4 categories there are 3 threshold parameters expected as a matrix --> rows: thresholds, columns: items)
}
  \item{theta0}{
starting value for the person parameter (default = 0)
}
  \item{exac}{
accuracy criterion - at which difference between two consecutive estimations the algorithm should stop
}
  \item{osch}{
maximum number of iterations (default is arbitrarily set to 33)
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an object of class "PPp" which contains the following values:

\item{resp}{entered response vector}
\item{estimate}{person parameter estimate}
\item{iterations}{number of iterations}
\item{estproc}{estimate for every iteration step - the resulting vector is of length(1:osch) - usually the algorithm stops before reaching the iteration limit, so there are NA's remaining at the end}
\item{thresh}{slightly modified threshold matrix (new first row added with zeroes only)}
\item{SE}{standard error}.
}
\references{
Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.

Muraki, Eiji (1992). A Generalized Partial Credit Model: Application of an EM Algorithm. Applied Psychological Measurement, 16, 159-176.

Muraki, Eiji (1993). Information Functions of the Generalized Partial Credit Model. Applied Psychological Measurement, 17, 351-363.

Wang, S. and Wang, T. (2001). Precision of Warm's Weighted Likelihood Estimates for a Polytomous Model in Computerized Adaptive Testing. Applied Psychological Measurement, 25, 317-331.
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
THRES  <- matrix(c(-2.41,0.36,-0.97,0.92,-0.65,-0.65,-0.28,0.72),nrow=2)
sl     <- c(0.5,1,1.5,1.1)
v1     <- c(1,2,1,3)

e1     <- PP_GPCMmle(u=v1,sl=sl,s=THRES)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
