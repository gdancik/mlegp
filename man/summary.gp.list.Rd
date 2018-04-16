\name{summary.gp.list}
\alias{summary.gp.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process List Summary Information }
\description{
  prints a summary of a Gaussian process list object, or (a subset) of its components
}
\usage{
\method{summary}{gp.list}(object, nums = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of type \code{gp.list} }
  \item{nums}{ optionally, a vector of integers corresponding to Gaussian processes in the list to summarize  }
  \item{\dots}{ for compatibility with generic method \code{summary} }
}
\details{
  if \code{nums} is \code{NULL}, prints out a summary of the Gaussian process list
 
  if \code{nums} is not \code{NULL}, displays a summary of the Gaussian processes specified by \code{nums} by calling \code{summary.gp} for each Gaussian process
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{summary.gp}} }
\examples{
x = -5:5 
y1 = sin(x) + rnorm(length(x),sd=.1)
y2 = sin(x) + 2*x + rnorm(length(x), sd = .1)
fitMulti = mlegp(x, cbind(y1,y2))
summary(fitMulti)  ## summary of the Gaussian process list
summary(fitMulti, nums = 1:2) ## summary of Gaussian processes 1 and 2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
\keyword{ methods }
