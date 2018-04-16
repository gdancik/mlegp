\name{mlegp-svd-functions}
\alias{mlegp-svd-functions}
\alias{pcweights}
\alias{getSingularValues}
\alias{numSingularValues}
\alias{singularValueImportance}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Singular Value Decomposition functions for mlegp }
\description{
	Functions that deal with the singular value decomposition of an output Y, for use with Gaussian process lists
}
\usage{
pcweights(Y, weights.num = NULL, cutoff = 99)
getSingularValues(Y)
singularValueImportance(Y)
numSingularValues(Y, cutoff = 99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ the output to decompose, where each column of \code{Y} correspond to k-dimensional output }
  \item{weights.num}{ optionally, the number of principle component weights to keep }
  \item{cutoff}{ if specified, \code{weights.num} is set to correspond to number of principle component weights that account for \code{cutoff} percent of the total variation in \code{Y}; if both \code{cutoff} and \code{weights.num} are specified, \code{weights.num} will be used}}
\details{
	
	Utilizes the singular value decomposition (SVD) of \code{Y}, Y = UDVprime. Columns of \code{Y} should correspond to a single k-dimensional observation (e.g., functional output of a computer model, evaluated at a particular input).      

	For a k x m matrix \code{Y}, and r = min(k,m), in the complete SVD, \code{U} is k x r, \code{D} is r x r, containing the singular values along the diagonal, and \code{Vprime} is r x m. The output \code{Y} is approximated by keeping l < r singular values, keeping a UD matrix of dimension k x l, and the \code{Vprime} matrix of dimension l x m. Each column of \code{Vprime} now contains l principle component weights, which can be used to reconstruct the functional output.        
}
\value{

  \code{pcweights} returns a list with components:
  \item{UD}{the UD matrix corresponding to the number of principle components kept}
  \item{Vprime}{The Vprime matrix corresponding to the number of principle components kept}

  Note: the number of principle component weights kept is equal to dim(UD)[2]

  \code{getSingularValues} returns a matrix containing the singular values of \code{Y}

  \code{numSingularValues} returns the minimum number of singular values accounting for \code{cutoff} percent of the variation in \code{Y}

  \code{singularValueImportance} returns a matrix where element i corresponds to the percentage of total variation in \code{Y} accounted for by the first i singular values 

}
\references{ Heitmann, K., Higdon, D., Nakhleh, C., Habib, S., 2006. Cosmic Calibration. The Astrophysical Journal, 646, 2, L1-L4.

\url{http://www1.easternct.edu/dancikg/} 
}

\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }

\note{ these functions are utilized by \code{mlegp} to fit Gaussian processes to principle component weights}
\seealso{ \code{\link{mlegp}}} 
\examples{


## create functional output that varies based on parameter 'p' ##
x = seq(-5,5,by=.2)
p = 1:50
y = matrix(0,length(p), length(x))
for (i in p) {
	y[i,] = sin(x) + i + rnorm(length(x), sd  = .1)
}

singularValueImportance(t(y))
numSingularValues(t(y), cutoff = 99.99)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ ts }
\keyword{ multivariate }

