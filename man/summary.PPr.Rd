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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
e1 <- PP_RMmle(c(-4,-3,-1.789,0.1,0.2,1.11,2.13),expol=TRUE)
summary(e1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
