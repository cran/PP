\name{PP-package}
\alias{PP-package}
\alias{PP}
\docType{package}
\title{
Estimation of Person Parameters for the 3-PL model and the GPCM under a ML and WL approach.
}
\description{
PP-package has been developed to easily compute ML and WL estimates (Warm 1989) of Person-Parameters for a given response vector and given Item Parameters of the 3-PL model (Birnbaum 1968) and the GPCM (Muraki 1992).
}
\details{
\tabular{ll}{
Package: \tab PP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-13\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Manuel Reif
}
\references{

Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In Lord, F.M. & Novick, M.R. (Eds.), Statistical theories of mental test scores. Reading, MA: Addison-Wesley.

Muraki, Eiji (1992). A Generalized Partial Credit Model: Application of an EM Algorithm. Applied Psychological Measurement, 16, 159-176.

Warm, Thomas A. (1989). Weighted Likelihood Estimation Of Ability In Item Response Theory. Psychometrika, 54, 427-450.
}

\keyword{ package }
\seealso{

}
\examples{
PP_3PLwle(u=c(1,0,1,0),a=c(1.1,2,0.2,1),s=c(-2,-1.23,1.11,3.48),i=c(0,0.1,0.3,0.2))
}
