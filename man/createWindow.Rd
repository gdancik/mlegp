\name{createWindow}
\alias{createWindow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process Plotting Functions }
\description{
  sets the number of sections on the current device based on the number of figures to draw
}
\usage{
createWindow(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ the number of figures to draw }
}
\details{

  Sets the graphical device so that the number of columns (ncol) is trunc(sqrt(n)) and the number of rows is ceiling(n/ncol)

This function is called by \code{mlegp} plotting functions that construct separate graphs for multiple parameters

}
\references{ \url{https://github.com/gdancik/mlegp/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\seealso{ This function may be used by \code{\link{plot.gp}}, among others }
\keyword{ internal }
\keyword{ device }

