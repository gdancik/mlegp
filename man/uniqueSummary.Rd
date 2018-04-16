\name{uniqueSummary}
\alias{uniqueSummary}
\title{ Summary of outputs for each unique input }
\description{
	Finds sample means and variances for a matrix of observations when replicates are present
}
\usage{
uniqueSummary(X,Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ matrix of inputs }
  \item{Y}{ matrix of outputs corresponding to x }
}
\details{
\code{uniqueSummary} calculates sample means and variances for each unique input. For input values with no replicates, the sample variance will be `NA'. This function is used by \code{mlegp} to fit a GP to a collection of means observations at each design point.  
}
\value{
	A list with the following components:
	\item{reps}{the number of reps for each unique design point}
	\item{uniqueX}{the input matrix with duplicate inputs removed}
	\item{uniqueMeans}{sample means corresponding to each unique input}
	\item{uniqueVar}{sample variances corresponding to each unique input}
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\examples{
x = c(1,1,2,2,3)
y = x + rnorm(length(x))
uniqueSummary(x,y)

}

\keyword{ methods }
