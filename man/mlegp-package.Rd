\name{mlegp-package}
\alias{mlegp-package}
\docType{package}
\title{
mlegp package
}
\description{
Maximum likelihood Gaussian process modeling for univariate and multi-dimensional outputs with diagnostic plots and sensitivity analysis. 
}
\details{
\tabular{ll}{
Package: \tab mlegp\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2007-12-05\cr
License: \tab Gnu General Public License (Version 3)\cr
}

This package obtains maximum likelihood estimates of Gaussian processes (GPs) for univariate and multi-dimensional outputs, for Gaussian processes with product exponential correlation structure; a constant or linear regression mean function; no nugget term, constant nugget term, or a nugget matrix that can be specified up to a multiplicative constant. The latter provides some flexibility for using GPs to model heteroscedastic responses.

Multi-dimensional output can be modelled by fitting independent GPs to each output, or to the most important principle component weights following singular value decomposition of the output. Plotting of main effects for functional output is also implemented.

Contact the maintainer for a package version that implements sensitivity analysis including Functional Analysis of Variance (FANOVA) decomposition, plotting functions to obtain diagnostic plots, main effects, and two-way factor interactions.

For a complete list of functions, use 'library(help="mlegp")'.

}
\author{ Garrett M. Dancik \email{ dancikg@easternct.edu } }
\references{

Santner, T.J. Williams, B.J., Notz, W., 2003. The Design and Analysis of Computer Experiments (New York: Springer).

Schonlau, M. and Welch, W. 2006. Screening the Input Variables to a Computer Model Via Analysis of Variance and Visualization, in Screening: Methods for Experimentation in Industry, Drug Discovery, and Genetics. A Dean and S. Lewis, eds. (New York: Springer). 
 
Heitmann, K., Higdon, D., Nakhleh, C., Habib, S., 2006. Cosmic Calibration. The Astrophysical Journal, 646, 2, L1-L4.

\url{https://github.com/gdancik/mlegp/} 

}
\keyword{ package }
