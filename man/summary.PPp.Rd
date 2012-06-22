\name{summary.PPp}
\alias{summary.PPp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes objects of class PPp
}
\description{
The results of the related functions are displayed in a reasonable way.
}
\usage{
\method{summary}{PPp}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class PPp
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

\author{
Manuel Reif
}

\examples{
THRES  <- matrix(c(-2.41,0.36,-0.97,0.92,-0.65,-0.65,-0.28,0.72),nrow=2)
sl     <- c(0.5,1,1.5,1.1)
v1     <- c(1,2,1,3)

e1     <- PP_GPCMmle(u=v1,sl=sl,s=THRES)
summary(e1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
