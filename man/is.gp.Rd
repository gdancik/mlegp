\name{is.gp}
\alias{is.gp}
\alias{is.gp.list}
\title{ Gaussian Process and Gaussian Process Lists }
\description{
     Test for a Gaussian process object or a Gaussian process list object
}
\usage{
is.gp(x)
is.gp.list(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object to be tested }
}
\value{
  \code{is.gp} returns \code{TRUE} or \code{FALSE}, depending on whether its argument inherits the \code{gp} class or not

  \code{is.gp.list} returns \code{TRUE} or \code{FALSE}, depending on whether its argument inherits the \code{gp.list} class or not
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\examples{
## fit a single Gaussian process ##
x = -5:5; y1 = sin(x) + rnorm(length(x),sd=.1)
fit1 = mlegp(x, y1)

is.gp(fit1)       ## returns TRUE
is.gp.list(fit1)  ## returns FALSE
}

\keyword{ methods }
