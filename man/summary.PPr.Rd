\name{summary.PPr}
\alias{summary.PPr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes objects of class PPr
}
\description{
The results of the related functions are displayed in a reasonable way.
}
\usage{
\method{summary}{PPr}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class PPr
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}


\author{
Manuel Reif
}


\examples{
e1 <- PP_RMmle(c(-4,-3,-1.789,0.1,0.2,1.11,2.13),expol=TRUE)
summary(e1)
}
