\name{plotObservedEffects}
\alias{plotObservedEffects}
\alias{plotObservedEffectsDefault}
\alias{plotObservedEffects.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Observed Values Vs. Each Dimension of the Design Matrix }
\description{
	Constructs multiple graphs, plotting each parameter from the design matrix on the x-axis and observations on the y-axis
}
\usage{
plotObservedEffects(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{gp} or a design matrix }
  \item{\dots}{ if x is a design matrix, a vector of observations;
		if x is of class \code{gp}, a vector of parameter numbers or parameter names to plot (by default, all parameters will be graphed)
		 }
}
\details{
	if \code{x} is NOT of class \code{gp} (i.e., \code{x} is a design matrix), all columns of \code{x} will be plotted separately against the vector of observations

	if \code{x} is of class \code{gp}, the specified columns of the design matrix of \code{x} will be plotted against the the observations
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{ It is often useful to use this function before fitting the gaussian process, to check that the observations
	are valid 
}
\examples{

## create the design and output matrices ##
x1 = kronecker(seq(0,1,by=.25), rep(1,5))
x2 = rep(seq(0,1,by=.25),5)
z = 4 * x1 - 2*x2 + x1 * x2 + rnorm(length(x1), sd = 0.001)

## look at the observed effects prior to fitting the GP ##
plotObservedEffects(cbind(x1,x2), z)

## fit the Gaussian process ##
fit = mlegp(cbind(x1,x2), z, param.names = c("x1", "x2"))

## look at the observed effects of the fitted GP (which are same as above)
plotObservedEffects(fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ methods }

