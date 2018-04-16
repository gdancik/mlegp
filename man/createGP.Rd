\name{createGP}
\alias{createGP}
\title{ creates a Gaussian process object }
\description{
	creates a Gaussian process \code{gp} object
}
\usage{
createGP(X, Z, beta, a, meanReg, sig2, nugget, 
         param.names = 1:dim(X)[2], constantMean = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ the design matrix }
  \item{Z}{ output obtained from the design matrix \code{X}, as a vector or a 1-column matrix }
  \item{beta}{ vector of correlation coefficients }
  \item{a}{ vector of smoothness parameters in the correlation function (if \code{mlegp} is used, these will be 2) }
  \item{meanReg}{ the constant mean if \code{constantMean} = 1, otherwise the regression coefficients of the mean function such that \code{meanReg} pre-multiplied by (1 X) will produce the mean matrix  }
  \item{sig2}{ the unconditional variance of the Gaussian process }
  \item{nugget}{ the constant nugget or a vector of length \code{nrow(X)} corresponding to the diagonal nugget matrix }
  \item{param.names}{ optional vector of parameter names (with length equal to \code{ncol(X)}}
  \item{constantMean}{ 1 if the Gaussian process has a constant mean; 0 otherwise }
}
\value{

   an object of class \code{gp} that contains the following components:

  \item{Z}{matrix of observations}
  \item{numObs}{ number of observations}
  \item{X}{the design matrix}
  \item{numDim}{ number of dimensions of X} 
  \item{constantMean}{1 if GP has a constant mean; 0 otherwise}
  \item{mu}{the mean matrix}
  \item{Bhat}{mean function regression coefficients}
  \item{beta}{correlation parameters}
  \item{a}{smoothness parameters in correlation function}
  \item{sig2}{unconditional variance of predicted expected output}
  \item{params}{vector of parameter names, corresponding to columns of \code{X}}
  \item{invVarMatrix}{inverse var-cov matrix}
  \item{nugget}{constant nugget or vector corresponding to the diagonal nugget matrix for a single observation generated from each element in \code{X} }
  \item{loglike}{the log likelihood of the observations}
  \item{cv}{results from cross-validation, where
	\code{cv[,1]} are the cross-validated predictions
  	\code{cv[,2]} are the variances of the cross-validated predictions}
}

\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{ this function is called by \code{mlegp} and should not be called by the user }
\seealso{ \code{\link{mlegp}} }
\keyword{ internal }
\keyword{ models }
\keyword{ smooth }
\keyword{ ts }
