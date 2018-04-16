\name{mlegp-parameter-lookup}
\alias{paramLookup}
\alias{toParamIndexes}
\alias{matchIndexes}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parameter Lookup Functions }
\description{
  These functions are used to match the name of a parameter with its position in a parameter list
}
\usage{
toParamIndexes(m, string)
matchIndexes(m, string)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{ a vector of names of the parameters of interest }
  \item{string}{ a vector of parameter names }
}
\value{
a vector where element i contains the position of \code{m[i]} in \code{string}

if \code{m} contains integers and \code{toParamIndexes} is called, \code{m} will be returned, without a check of whether
or not the indices are valid

if \code{m[i]} is not an element of \code{string}, \code{toParamIndexes} will display an error, whereas \code{matchIndexes} will return \code{NA}

}
\references{ \url{http://www1.easternct.edu/dancikg/} }
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{  
  this function does not need to be called by the user; it exists so that the user can pass in a vector of parameter numbers or parameter names to various functions when specifying a subset of the parameters of the Gaussian process design matrix
  
for \code{toParamIndexes}, \code{m} can contain integers or characters, but cannot contain both 
}
\examples{

param.names = c("one", "two", "three")
toParamIndexes(c("one", "three"), param.names)
#toParamIndexes(c("four"), param.names) # will give an error
}

\keyword{ internal }
