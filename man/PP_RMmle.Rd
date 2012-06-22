\name{PP_RMmle}
\alias{PP_RMmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of ML Person Parameters for the Rasch Model
}
\description{
This is a wrapper function, which calls PP_3PLmle to estimate person parameters for items calibrated accordingly to the Rasch Model. So just the item difficulty parameters have to be declared, the slope is assumed to be 1 for all items - the lower-asymtote is set to 0. The main result is a vector which contains all possible person parametes for a given set of items. If required, the algorithm extrapolates for the marginal totals of zero, and all items solved.
}

\usage{
PP_RMmle(s, expol = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
numeric vector of difficulty parameters (item difficulty)
}
  \item{expol}{
TRUE if extrapolation for minimum and maximum raw scores (0 and all items solved) should be performed, FALSE otherwise. If the number of items is less then 4, a standard linear model (OLS) is applied to perform extrapolation - otherwise a cubic spline (package: splines) is used.
}
  \item{\dots}{
refers to further arguments of PP_3PLmle
}
}


\value{
Returns an object of class "PPr" which contains the following values:

\item{result}{contains the whole list of output returned from PP_3PLmle, with all the known details for each raw-score}
\item{estimates}{matrix with raw score, person parameter estimates and standard error }
\item{extrap}{returns whether and which kind of extrapolation was performed}

}
\references{
Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.

}
\author{
Manuel Reif
}



\seealso{
\code{\link{PP_3PLmle}},
\code{\link{PP_3PLwle}}
}

\examples{
PP_RMmle(c(-4,-3,-1.789,0.1,0.2,1.11,2.13),expol=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
