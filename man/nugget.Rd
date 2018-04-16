\name{mlegp-nugget-related}
\alias{nugget}
\alias{varPerReps}
\alias{estimateNugget}
\alias{anyReps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian Process Nugget Related Functions }
\description{
  Functions for detecting replicates and for calculating sample variance at specific design points
}
\usage{
varPerReps(X, Y)
estimateNugget(X, Y)
anyReps(X)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ the design matrix }
  \item{Y}{ a vector (or 1 column matrix) of observations }
}
\value{
  \code{varPerReps} returns a 1-column matrix where element i corresponds to the sample variance in observations corresponding to design point X[i]


  \code{estimateNugget} returns a double calculated by taking the mean of the matrix returned by \code{varPerReps}


  \code{anyReps} returns \code{TRUE} if two or more rows of \code{X} are identical
}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{  
These functions are used by \code{mlegp} to set an initial value of the nugget when a constant nugget is being estimated. The function \code{varPerReps} may also be useful for specifying the form of the nugget matrix for use with \code{mlegp}. 
}
\examples{

x = matrix(c(1,1,2,3,3))   # the design matrix
y = matrix(c(5,6,7,0,10))  # output

anyReps(x)
varPerReps(x,y)
estimateNugget(x,y)

}

\keyword{ htest }
