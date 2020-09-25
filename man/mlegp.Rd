\name{mlegp}
\alias{mlegp}
\alias{mlegp2}
\alias{fitGPfromR}
\title{ mlegp: maximum likelihood estimation of Gaussian process parameters }
\description{
  Finds maximum likelihood estimates of Gaussian process parameters for a vector (or matrix) of one (or more) responses. For multiple responses, the user chooses between fitting independent Gaussian processes to the separate responses or fitting independent Gaussian processes to principle component weights obtained through singular value decomposition of the output. The latter is useful for functional output or data rich situations. 
}
\usage{
mlegp(X, Z, constantMean = 1, nugget = NULL, nugget.known = 0, 
      min.nugget = 0, param.names = NULL, gp.names = NULL, 
	  PC.UD = NULL, PC.num = NULL, PC.percent = NULL, 
	  simplex.ntries = 5, simplex.maxiter = 500, simplex.reltol = 1e-8,  
	  BFGS.maxiter = 500, BFGS.tol = 0.01, BFGS.h = 1e-10, seed = 0, 
      verbose = 1, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ the design matrix }
  \item{Z}{ vector or matrix of observations; corresponding to the rows of \code{X} }
  \item{constantMean}{ a value of 1 indicates that each Gaussian process will have a constant mean; otherwise the mean function will be a linear regression in \code{X}, plus an intercept term }
  \item{nugget}{ if nugget.known is 1, a fixed value to use for the nugget or a vector corresponding to the fixed diagonal nugget matrix; otherwise, either a positive initial value for the nugget term which will be estimated, or a vector corresponding to the diagonal nugget matrix up to a multiplicative constant. If \code{NULL} (the default), \code{mlegp} estimates a nugget term only if there are replicates in the design matrix, see details}
  \item{nugget.known}{1 if a plug-in estimate of the nugget will be used; 0 otherwise}
  \item{min.nugget}{ minimum value of the nugget term; 0 by default }
  \item{param.names}{ a vector of parameter names, corresponding to the columns of X; parameter names are
			\sQuote{p1}, \sQuote{p2}, ... by default }
  \item{gp.names}{ a vector of GP names, corresponding to the GPs fit to each column of \code{Z} or each PC weight }
  \item{PC.UD}{ the UD matrix if \code{Z} is a matrix of principle component weights; see \code{\link{mlegp-svd-functions}} }
  \item{PC.num}{ the number of principle component weights to keep in the singular value decomposition of \code{Z} }
  \item{PC.percent}{ if not \code{NULL} the number of principle component weights kept is the minimum number that accounts for \code{PC.percent} of the total variance of the matrix \code{Z} }
  \item{simplex.ntries}{ the number of simplexes to run }
  \item{simplex.maxiter}{ maximum number of evaluations / simplex }
  \item{simplex.reltol}{ relative tolerance for simplex method, defaulting to 1e-16 }
  \item{BFGS.maxiter}{ maximum number of iterations for BFGS method }
  \item{BFGS.tol}{ stopping condition for BFGS method is when norm(gradient) < BFGS.tol * max(1, norm(x)), where x is the parameter vector and norm is the Euclidian norm }
  \item{BFGS.h}{ derivatives are approximated as [f(x+BFGS.h) - f(x)] / BFGS.h)  }
  \item{seed}{ the random number seed }
  \item{verbose}{ a value of '1' or '2' will result in status updates being printed; a value of '2' results in more information }
  \item{parallel}{ if TRUE will fit GPs in parallel to each column of Z, or each set of PC weights; See details }
}
\details{
This function calls the C function fitGPFromR which in turn calls fitGP (both in the file fit_gp.h) to fit each Gaussian process.


Separate Gaussian processes are fit to the observations in each column of \code{Z}. Maximum likelihood estimates for correlation and nugget parameters are found through numerical methods (i.e., the Nelder-Mead Simplex and the L-BFGS method), while maximum likelihood estimates of the mean regression parameters and overall variance are calculated in closed form (given the correlation and (scaled) nugget parameters). Multiple simplexes are run, and estimates from the best simplex are used as initial values to the gradient (L-BFGS) method. 

Gaussian processes are fit to principle component weights by utilizing the singular value decomposition (SVD) of \code{Z}, Z = UDVprime. Columns of \code{Z} should correspond to a single k-dimensional observation (e.g., functional output of a computer model, evaluated at a particular input)      

In the complete SVD, \code{Z} is k x m, and r = min(k,m), \code{U} is k x r, \code{D} is r x r, containing the singular values along the diagonal, and \code{Vprime} is r x m. The output \code{Z} is approximated by keeping l < r singular values, keeping a UD matrix of dimension k x l, and the \code{Vprime} matrix of dimension l x m. Each column of \code{Vprime} now contains l principle component weights, which can be used to reconstruct the functional output.        

If nugget.known = 1, nugget = NULL, and replicate observations are present, the nugget will be fixed at its best linear unbiased estimate (a weighted average of sample variances). 
For each column of Z, a GP will be fit to a collection of sample means rather than all observations. This is the recommended approach as it is more accurate and computationally more efficient. 

Parallel support is provided through the package \code{snowfall} which allows multiple GPs to be fit in parallel. The user must set up the cluster using \code{sfInit} and call \code{sfLibrary(mlegp)} to load the library onto the slave nodes. Note: GP fitting is not recommended when the number of observations are large (> 100), in which case sequential GP fitting is faster. } 

\value{
  an object of class \code{gp.list} if Z has more than 1 column, otherwise an object of class \code{gp} 
}
\references{ 

Santner, T.J. Williams, B.J., Notz, W., 2003. The Design and Analysis of Computer Experiments (New York: Springer).

Heitmann, K., Higdon, D., Nakhleh, C., Habib, S., 2006. Cosmic Calibration. The Astrophysical Journal, 646, 2, L1-L4.

Dancik, GM and Dorman, KS (2008). mlegp: statistical analysis for
computer models of biological systems using R. Bioinformatics 24(17), pp. 1966-1967

\url{https://github.com/gdancik/mlegp/} 

}
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\note{ The random number seed is 0 by default, but should be randomly set by the user

In some situations, especially for noiseless data, it may be desirable to force a nugget term in order to make the variance-covariance matrix of the Gaussian process more stable; this can be done by setting the argument \code{min.nugget}.

If fitting multiple Gaussian processes, the arguments \code{min.nugget} and \code{nugget} apply to all Gaussian processes being fit.

In some cases, the variance-covariance matrix is stable in C but not stable in R. When this happens, this function will attempt to impose a minimum value for the nugget term, and this will be reported. However, the user is encouraged to refit the GP and manually setting the argument \code{min.nugget} in \code{mlegp}.

When fitting Gaussian processes to principle component weights, a minimum of two principle component weights must be used.

}
\seealso{ \code{\link{createGP}} for details of the \code{gp} object; \code{\link{gp.list}} for details of the \code{gp.list} object; \code{\link{mlegp-svd-functions}} for details on fitting Gaussian processes to high-dimensional data using principle component weights; the L-BFGS method uses C code written by Naoaki Okazaki (http://www.chokkan.org/software/liblbfgs)  

}
\examples{

###### fit a single Gaussian process ######
x = -5:5; y1 = sin(x) + rnorm(length(x),sd=.1)
fit1 = mlegp(x, y1)

## summary and diagnostic plots ##
summary(fit1)
plot(fit1)

###### fit a single Gaussian process when replciates are present ######
x = kronecker(-5:5, rep(1,3))
y = x + rnorm(length(x))

## recommended approach: GP fit to sample means; nugget calcualted from sample variances ##
fit1 = mlegp(x,y, nugget.known = 1)

## original approach: GP fit to all observations; look for MLE of nugget ##
fit2 = mlegp(x,y)


###### fit multiple Gaussian processes to multiple observations ######
x = -5:5 
y1 = sin(x) + rnorm(length(x),sd=.1)
y2 = sin(x) + 2*x + rnorm(length(x), sd = .1)
fitMulti = mlegp(x, cbind(y1,y2))

## summary and diagnostic plots ##
summary(fitMulti)
plot(fitMulti)


###### fit multiple Gaussian processes using principle component weights ######

## generate functional output ##
x = seq(-5,5,by=.2)
p = 1:50
y = matrix(0,length(p), length(x))
for (i in p) {
	y[i,] = sin(x) + i + rnorm(length(x), sd  = .01)
}

## we now have 10 functional observations (each of length 100) ##
for (i in p) {
	plot(x,y[i,], type = "l", col = i, ylim = c(min(y), max(y)))
	par(new=TRUE)
}

## fit GPs to the two most important principle component weights ##
numPCs = 2
fitPC = mlegp(p, t(y), PC.num = numPCs)
plot(fitPC) ## diagnostics

## reconstruct the output Y = UDV'
Vprime = matrix(0,numPCs,length(p))
Vprime[1,] = predict(fitPC[[1]])
Vprime[2,] = predict(fitPC[[2]])

predY = fitPC$UD\%*\%Vprime
m1 = min(y[39,], predY[,39])
m2 = max(y[39,], predY[,39])

plot(x, y[39,], type="l", lty = 1, ylim = c(m1,m2), ylab = "original y" )
par(new=TRUE)
plot(x, predY[,39], type = "p", col = "red", ylim = c(m1,m2), ylab = "predicted y" )

\dontrun{
### fit GPs in parallel ###
library(snowfall)
sfInit(parallel = TRUE, cpus = 2, slaveOutfile = "slave.out")
sfLibrary(mlegp)
fitPC = mlegp(p, t(y), PC.num = 2, parallel = TRUE)
}

}

\keyword{ models }
\keyword{ smooth }
\keyword{ ts }
