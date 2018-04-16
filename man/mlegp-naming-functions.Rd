\name{mlegp-naming-functions}
\alias{mlegp-naming-functions}
\alias{setParams}
\alias{setGPNames}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ mlegp naming functions }
\description{
	Functions that set the design matrix parameter names in a Gaussian process or Gaussian process list or the names of the Gaussian processes in a Gaussian process list
}
\usage{
setParams(x, s)
setGPNames(x, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of type \code{gp} or \code{gp.list} }
  \item{s}{ a vector of parameter names or Gaussian process names to set }
}
\details{

	\code{setParams} sets the parameter names of a Gaussian process (\code{gp}) object or of all Gaussian processes in a Gaussian process list (\code{gp.list}) object. \code{setGPNames} sets the names of the Gaussian processes in an object of type \code{gp.list}
}
\value{
  the object \code{x} with parameter or Gaussian process names set   
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{
Parameter and Gaussian process names are used in the output of various plotting functions. Both of these can also be set when the Gaussian process (list) is created by \code{\link{mlegp}}
}
\seealso{ \code{\link{mlegp}} for setting parameter and Gaussian process names during object creation }

\examples{
## fit multiple Gaussian processes to multiple observations ##
x = -5:5
y1 = sin(x) + rnorm(length(x),sd=.1)
y2 = sin(x) + 2*x + rnorm(length(x), sd = .1)
fitMulti = mlegp(x, cbind(y1,y2))

## plot diagnostics with default gp names ##
plot(fitMulti)

## change names and plot again ##
fitMulti = setGPNames(fitMulti, c("y1", "y2"))
plot(fitMulti)

## plot diagnostic for the first Gaussian process, predicted vs. parameter ##
plot(fitMulti[[1]], type = 2)

## change parameter names (of all Gaussian processes) and plot again ##
fitMulti = setParams(fitMulti, "param 1")
plot(fitMulti[[1]], type = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

