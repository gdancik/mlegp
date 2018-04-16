\name{CV}
\alias{CV}
\alias{calcPredictionErrorCV}
\alias{predictNewYCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian process cross-validation  }
\description{
  For a Gaussian process, calculates cross-validated predictions and the variance of cross-validated predictions for all points of the design. These are cross-validated in the sense that when predicting output at design point x, all observations at x are removed from the collection of observed outputs
}
\usage{
CV(gp, predictObserved = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gp}{ an object of type \code{gp} }
  \item{predictObserved}{ if \code{TRUE} (the default), the variance of a predicted observation is calculated (and a nugget term is added);  if \code{FALSE}, the variance of a predicted mean response is calculated (and a nugget term is NOT added)
 }
  \item{verbose}{ if \code{TRUE}, status is displayed following the cross-validation of each observation }

}
\value{
  a matrix where the first column corresponds to the cross-validated predictions and the second column corresponds to the variance of the cross-validated predictions
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ \code{\link{predict.gp}} }
\examples{

## fit a single Gaussian process ##
x = -5:5; y1 = sin(x) + rnorm(length(x),sd=.1)
fit1 = mlegp(x, y1)

cv = CV(fit1)  ## note that cv is the same as fit1$cv

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal  }
\keyword{ models }
