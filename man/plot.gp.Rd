\name{plot.gp}
\alias{plot.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Diagnostic Plots for Gaussian processes }
\description{

  Cross-Validated Diagnostic Plots for Gaussian Processes
}
\usage{
\method{plot}{gp}(x, type = 0, params = NULL, sds = 1, CI.at.point = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{gp} }
  \item{type}{ the type of graph to plot, 0 by default (see Details) }
  \item{params}{ for graph types 2 and 3, a vector of parameter names (or parameter indices) to plot against. By default, all parameters are looked at }
  \item{sds}{ the number of standard deviations to use for confidence bands/intervals, for graph types 0-3 }
  \item{CI.at.point}{ if TRUE, will plot confidence intervals at each predicted point, rather than bands, which is the default}
  \item{\dots}{ additional arguments to plot, but cannot overwrite xlab or ylab}
}
\details{
All plots involve cross-validated predictions and/or cross-validated standardized residuals. The cross-validation is in the sense that for predictions made at design point x, all observations at design point x are removed from the training set. 

Where relevant, open circles correspond to Gaussian process predictions, black lines correspond to the observations, and red lines correspond to confidence bands. The argument \code{type} determines the type of graph displayed, and is one of the following integers:\cr

\tabular{ll}{
		\tab 0 for observed vs. predicted AND observed vs. standardized residual (default),\cr 
		\tab 1 for observed vs. predicted only,\cr
		\tab 2 for parameter vs. predicted for all parameters,\cr
		\tab 3 for parameter vs. standardized residual for all parameters,\cr
		\tab 4 for normal quantile plot and histogram of standardized residuals } 
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{CV}} for cross-validation, \code{\link{plot.gp.list}} for plotting \code{gp.list} objects }
\examples{

## fit the gp ##
x = seq(-5,5,by=.5)
y = sin(x) + rnorm(length(x), sd=.1)
fit = mlegp(x,y)

## plot diagnostics ##
plot(fit)
plot(fit, type = 2)
}

\keyword{ hplot }
\keyword{ methods }
