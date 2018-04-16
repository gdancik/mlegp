\name{gp.list}
\alias{gp.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process Lists }
\description{
	Creates an object of type \code{gp.list}, given a list of Gaussian processes fit to separate sets of observations, or a list of Gaussian processes fit to principle component weights to approximate output of high dimension
}
\usage{
gp.list(..., param.names = NULL, UD = NULL, gp.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ either a \code{list} object, where each element is a Gaussian process; or several objects of type \code{gp} }
  \item{param.names}{ optionally, the parameter names corresponding to the columns of the design matrix of all Gaussian processes. By default, this will be equal to the parameter names of the first Gaussian process in \code{\dots}}
  \item{UD}{ the UD matrix, if the Gaussian process is fit to principle component weights }
  \item{gp.names}{ optionally, a vector of names for the Gaussian processes, defaulting to \sQuote{gp \#1}, \sQuote{gp \#2}, ...
}
}

\value{
  
  A \code{gp.list} object is a \code{list} object, where the first k elements correspond to k Gaussian processes passed in as \code{\dots}. This makes it straightforward to access a single Gaussian process. In addition, \code{gp.list} contains components:
  \item{params}{a vector of parameter names, corresponding to the columns of the design matrix}
  \item{numGPs}{the number of Gaussian processes in the list}
  \item{numDim}{the number of parameters in the design matrix}
  \item{numObs}{the number of observations}
  \item{names}{the names of the Gaussian processes}
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{ 
currently, we require that all Gaussian processes have the same dimension (number of columns in the design matrix) and the same number of observations

this function is called by \code{mlegp} and should not be called by the user
 }
\seealso{ \code{\link{mlegp}}, \code{\link{mlegp-svd-functions}} for more details about the UD matrix }
\examples{
  x = -5:5
  y1 = sin(x) + rnorm(length(x), sd=.1)
  y2 = sin(x) + rnorm(length(x), sd = .5)

  ## create the gp.list object ## 
  fitMulti = mlegp(x, cbind(y1,y2))

  plot(fitMulti)
 
  fitMulti   ## print summary of of the fitted Gaussian process list
  fitMulti[[2]]  ## print summary for the 2nd Gaussian process 

}

\keyword{ internal }
\keyword{ list }
\keyword{ models }
\keyword{ smooth }
\keyword{ ts }
