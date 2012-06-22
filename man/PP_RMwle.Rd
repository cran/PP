\name{PP_RMwle}
\alias{PP_RMwle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of WL Person Parameters for the Rasch Model
}
\description{
This is a wrapper function, which calls PP_3PLwle to estimate person parameters for items calibrated accordingly to the Rasch Model under Warms weighted maximum-likelihood approach (Warm 1989). So just the item difficulty parameters have to be declared, the slope is assumed to be 1 for all items - the lower-asymtote is set to 0. The main result is a matrix containing the estimated person parameters corresponding to all possible raw scores for a given set of items.
}
\usage{
PP_RMwle(s, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
numeric vector of difficulty parameters (item difficulty)
}
  \item{\dots}{
refers to further arguments of PP_3PLwle
}
}


\value{
Returns an object of class "PPr" which contains the following values:

\item{result}{contains the whole list of output returned from PP_3PLmle, with all the known details for each raw-score}
\item{estimates}{matrix with raw score, person parameter estimates and standard error }

}
\references{
Baker, Frank B., and Kim, Seock-Ho (2004). Item Response Theory - Parameter Estimation Techniques. CRC-Press.

Warm, Thomas A. (1989). Weighted Likelihood Estimation Of Ability In Item Response Theory. Psychometrika, 54, 427-450.
}
\author{
Manuel Reif
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PP_3PLwle}},
\code{\link{PP_3PLmle}}
}

\examples{
PP_RMwle(c(-4,-3,-1.789,0.1,0.2,1.11,2.13))
}
