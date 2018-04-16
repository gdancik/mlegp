\name{plot.gp.list}
\alias{plot.gp.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Diagnostics Plots for Gaussian Process Lists }
\description{

  Cross-validated Diagnostic Plots For Gaussian Process Lists

}
\usage{
\method{plot}{gp.list}(x, sds = 1, CI.at.point = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{gp.list} }
  \item{sds}{ the number of standard deviations to use for confidence bands / intervals}
  \item{CI.at.point}{ if TRUE, will plot confidence intervals around each predicted point, rather than bands, which is the default}
  \item{\dots}{ not used; for compatibility with \code{plot.gp} }
}
\details{
All plots involve cross-validated predictions and/or cross-validated standardized residuals. The cross-validation is in the sense that for predictions made at design point \code{x}, all observations at design point \code{x} are removed from the training set. 

Where relevant, open circles correspond to Gaussian process output, black lines correspond to the observations, and red lines correspond to confidence bands.

For each Gaussian process in \code{x}, \code{plot.gp} is called using graph type 1, which plots cross-validated predictions vs. observed values.
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{plot.gp}}, \code{\link{CV}} }
\examples{

## create data for multiple responses ##
x = seq(-5,5)
z1 = 10 - 5*x + rnorm(length(x))
z2 = 4 - 5*x + rnorm(length(x))
z3 = 7*sin(x) + rnorm(length(x))

## fit multiple Gaussian processes ##
fitMulti = mlegp(x, cbind(z1,z2,z3))

## plot diagnostics ##
plot(fitMulti)

}

\keyword{ hplot }
\keyword{ methods }
