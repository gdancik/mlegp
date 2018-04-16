\name{print.gp.list}
\alias{print.gp.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process List Summary Information }
\description{
  prints a summary of a Gaussian process list object, or (a subset) of its components
}
\usage{
\method{print}{gp.list}(x, nums = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{gp.list} }
  \item{nums}{ optionally, a vector of integers corresponding to Gaussian processes in the list to summarize  }
  \item{\dots}{ for compatibility with generic method \code{print} }
}
\details{
  if \code{nums} is \code{NULL}, prints summary information for the Gaussian process list object \code{x}, using \code{summary.gp.list}

  if \code{nums} is not \code{NULL}, prints summary information for each Gaussian process specified by \code{nums}, using \code{summary.gp}
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{summary.gp.list}}, \code{\link{summary.gp}} 
}
\examples{
x = -5:5 
y1 = sin(x) + rnorm(length(x),sd=.1)
y2 = sin(x) + 2*x + rnorm(length(x), sd = .1)
fitMulti = mlegp(x, cbind(y1,y2))
print(fitMulti)  ## summary of the Gaussian process list
print(fitMulti, nums = 1:2) ## summary of Gaussian processes 1 and 2

}
\keyword{ print }
\keyword{ methods }
