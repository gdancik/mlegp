\name{predict.gp}
\alias{predict.gp}
\alias{predictMu}
\alias{predictNewZ}
\alias{calcPredictionError}
\alias{predict.gp.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process Predictions }
\description{

  Gaussian Process Predictions
}
\usage{
\method{predict}{gp}(object, newData = object$X, se.fit = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an object of class \code{gp} }
  \item{newData}{ an optional data frame or matrix with rows corresponding to inputs for which to predict.
                  If omitted, the design matrix \code{X} of \code{object} is used. }
  \item{se.fit}{ a switch indicating if standard errors are desired }
  \item{\dots}{ for compatibility with generic method \code{predict} }
}
\details{
  
  The Gaussian process is used to predict output at the design points \code{newData}; if the logical \code{se.fit} is set to \code{TRUE}, standard errors (standard deviations of the predicted values) are also calculated. Note that if the Gaussian process contains a nugget term, these standard deviations correspond to standard deviations of predicted expected values, and NOT standard deviations of predicted observations. However, the latter can be obtained by noting that the variance of a predicted observation equals the variance of the predicted expected value plus a nugget term.

  If \code{newData} is equal to the design matrix of \code{object} (the default), and there is no nugget term, the Gaussian process interpolates the observations and the
  predictions will be identical to component \code{Z} of \code{object}. For cross-validation, the function \code{CV} should be used. 

}
\value{
  \code{predict.gp} produces a vector of predictions. If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
  \item{fit}{vector as above}
  \item{se.fit}{standard error of the predictions}
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{ for predictions with \code{gp.list} objects, call \code{predict.gp} separately for each gp in the list
}
\seealso{For cross-validated predictions, see \code{\link{CV}}  }
\examples{

x = -5:5; y = sin(x) + rnorm(length(x), sd = 0.001)
fit = mlegp(x,y)
predict(fit, matrix(c(2.4, 3.2)))
## predictions at design points match the observations 
## (because there is no nugget)
round(predict(fit) - fit$Z, 6)   

# this is not necessarily true if there is a nugget
fit = mlegp(x,y, min.nugget = 0.01)
round(predict(fit) - fit$Z,6)   

}

\keyword{ models }
\keyword{ smooth }
