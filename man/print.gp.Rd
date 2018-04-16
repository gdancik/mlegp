\name{print.gp}
\alias{print.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process Summary Information }
\description{
  prints a summary of a Gaussian process object
}
\usage{
\method{print}{gp}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{gp} }
  \item{\dots}{ for compatibility with generic method \code{print} }
}
\details{
  prints a summary of the Gaussian process object \code{x}, by calling \code{summary.gp}
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{summary.gp} for more description of the output} }
\examples{

x = -5:5; y1 = sin(x) + rnorm(length(x),sd=.1)
fit1 = mlegp(x, y1)
print(fit1)

}
\keyword{ print }
\keyword{ methods }


