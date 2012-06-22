\name{PP_GPCMmle}
\alias{PP_GPCMmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of ML Person Parameters for the GPCM
}
\description{
This function estimates Person Parameters for a given response vector (polychotom) and fixed parameters ("slope" and threshold-parameters) of the generalized partial credit model (Muraki 1992) under a maximum-likelihood approach.
}
\usage{
PP_GPCMmle(u, sl, s, theta0 = 0, exac = 0.001, osch = 33, ...)
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
Note that if a set of items has unequal number of categories, then "missing" categories (items with a smaller number of categories) have to be considered as missing values in the matrix of thresholds s (see example below).
}


\value{
Returns an object of class "PPp" which contains the following values:

\item{resp}{entered response vector}
\item{estimate}{person parameter estimate}
\item{iterations}{number of iterations}
\item{estproc}{estimate for every iteration step - the resulting vector is of length(1:osch) - usually the algorithm stops before reaching the iteration limit, so there are NA's remaining at the end}
\item{thresh}{slightly modified threshold matrix (new first row added with zeroes only)}
\item{SE}{standard error}

}
\references{

Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.

Muraki, Eiji (1992). A Generalized Partial Credit Model: Application of an EM Algorithm. Applied Psychological Measurement, 16, 159-176.

}
\author{
Manuel Reif
}


\seealso{
\code{\link{PP_GPCMwle}}
}

\examples{

THRES  <- matrix(c(-2.41,0.36,-0.97,0.92,-0.65,-0.65,-0.28,0.72),nrow=2)
sl     <- c(0.5,1,1.5,1.1)
v1     <- c(1,2,1,3)

e1     <- PP_GPCMmle(u=v1,sl=sl,s=THRES)


# an ordinary PCM with a different number of categories

THRES2  <- matrix(c(-2.41,NA,-0.57,0.92,-0.65,-0.65,-0.28,1.19),nrow=2)
sl2     <- c(1,1,1,1)
v2     <- c(1,2,1,3)

e2     <- PP_GPCMmle(u=v2,sl=sl2,s=THRES2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
