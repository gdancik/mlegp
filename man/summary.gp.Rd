\name{summary.gp}
\alias{summary.gp}
\title{ Gaussian Process Summary Information }
\description{
  prints a summary of a Gaussian process object 
}
\usage{
\method{summary}{gp}(object, ...)
}
\arguments{
  \item{object}{ an object of class \code{gp} }
  \item{\dots}{ for compatibility with generic method \code{summary} } 
}
\details{
  prints a summary of the Gaussian process object \code{object}. Output should be self explanatory, except for possibly CV RMSE, the cross-validated root mean squared error (the average squared difference between the observations and cross-validated predictions); and CV RMaxSE, the maximum cross-validated root squared error. If the design in the Gaussian process object contains any replicates, the root mean pure error (RMPE), which is the square root of the average within replicate variance and the root max pure error (RMaxPE) are also reported.   
  
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{createGP}} for details of the Gaussian process object
}
\examples{

## no replicates in the design matrix ##
x1 = -5:5; y1 = sin(x1) + rnorm(length(x1),sd=.1)
fit1 = mlegp(x1, y1)
summary(fit1)

## with replicates in the design matrix ##
x2 = kronecker(x1, rep(1,3)) 
y2 = sin(x2) + rnorm(length(x2), sd = .1)
fit2 = mlegp(x2,y2)
summary(fit2)

}
\keyword{ print }
\keyword{ methods }
