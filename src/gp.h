
/***************************************************************************
    Maximum Likelihood Estimation of Gaussian Processes (mlegp)
    Copyright (C) 2007, Garrett Dancik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    For questions or comments, contact Garrett Dancik (dancikg@nsula.edu)
****************************************************************************/

/*******************************************************************************
gp.h - finds mles of gaussian process parameters that have closed form
	solutions (i.e, B in the mean matrix X'B, and the overall variance sig2;
	also creates a correlation matrix assuming product exponential 
	structure
*******************************************************************************/

#ifndef __GP__
#define __GP__

#include <math.h>
#include <float.h>
#include "matrix_vector.h"

#ifndef M_LNPI
#define M_LNPI 1.1447298858494
#endif

#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif

/*********************************************************************************
calcBhat - finds the mle of the mean regression term X'B of the gaussian process
	double *X - the design matrix for the regression term X'B
	int nrowX - number of rows of X (same as # of observations and length 
			width of Vinv) 
	int ncolX - number of columns of X
	double *Vinv - the inverse variance-covariance matrix specified up to
			   a multiplicative constant
	const double *Y - the observations    
	double *bhat - vector that stores the mle of B (that has length ncolX)
	
	return value: 1 if there is an error inverting X'X; 0 otherwise

	notes: X is for the mean regression function, not the design matrix for 
		comp model runs
**********************************************************************************/
int calcBhat(const double *X, const int nrowX, const int ncolX, const double *Vinv,  
	     const double *Y, double *bhat) {

	double *Xprime = MATRIX(nrowX, ncolX);
	transpose(X, Xprime, nrowX, ncolX);

	double *XprimeVinv = MATRIX(ncolX, nrowX);   // ncolX is nrow of Xprime
	matrix_multiply(Xprime, ncolX, nrowX, Vinv, nrowX, XprimeVinv);

	double *XprimeVinvX = MATRIX(ncolX, ncolX);

	matrix_multiply(XprimeVinv, ncolX, nrowX, X, ncolX, XprimeVinvX);	

	double *ans1 = MATRIX(ncolX, ncolX);
	createIdentityMatrix(ans1, ncolX);

	int *piv = malloc(ncolX*sizeof (int));
	int info = LP_gen_solve(XprimeVinvX, ncolX, ans1, ncolX, piv);
	if (info != 0) {
		FREE_MATRIX(Xprime);
		FREE_MATRIX(XprimeVinv);
		FREE_MATRIX(ans1);
		free(piv);
		return 1;
	}

	double *ans2 = MATRIX(ncolX, 1);
	matrix_multiply(XprimeVinv,ncolX, nrowX, Y, 1, ans2);
	matrix_multiply(ans1, ncolX, ncolX, ans2, 1, bhat);

	free(piv);
	FREE_MATRIX(Xprime);
	FREE_MATRIX(XprimeVinv);
	FREE_MATRIX(XprimeVinvX);
	FREE_MATRIX(ans1);
	FREE_MATRIX(ans2);
	return 0;
}


/**************************************************************************************
calcMLEsig2 - finds the mle of the overall variance of the GP;
	const double *Y - the observations   
	const double *mu - mean matrix with same dimensions as Y 
	const double *Vinv - the inverse variance-covariance matrix to mult. constant
	return value: the mle of sig2
**************************************************************************************/

double calcMLESig2 (const double *Y, const double *mu, const double *Vinv, int n) {

	double *diff = VECTOR(n);
	vectorCopy(Y, diff, n);
	vectorSubtract(diff,mu, n);	// Y - mu

	double *diffVinv = MATRIX(1,n);
	xprimeA(diff, Vinv, diffVinv, n,n);
	
	double sig2 = dotprod(diffVinv, diff,n) / (double) n;
		
	FREE_VECTOR(diff);
	FREE_MATRIX(diffVinv);	
	return sig2;
}

// make sure corr is really correlations 
// assume that nreps / design point is the SAME, and nY is # of design points
// nY MUST be equal to # of design points for calculation to be correct
// returns an unbiased estimate of Sigma2(GP)
// requires that MSE < MSA
/***********************
double calcANOVASig2 (double MSA, double MSE, const double *corr, const double nY, const double nreps) { 
	double num = (MSA - MSE) * (nY - 1.00);
	double den = nreps*(nY + 1.0 - 2.0 / nY * sumPackedMatrix(corr, nY)); 
	return num / den;
}
****************************/

/******************************************************************************************
addNugget - adds a constant nugget term to the diagonal of the packed matrix m that is nxn 
******************************************************************************************/
void addNuggetToPackedMatrix(double *m, double nugget, int n) {
	int i = 0;
	while (n > 0) {
		m[i] += nugget;
		i+= n;
		n--;
	}
}


/*********************************************************************************************
addNuggetMatrixToPackedMatrix - adds c * diag(nugget matrix) to the packed matrix m that is nxn
**********************************************************************************************/
void addNuggetMatrixToPackedMatrix(double *m, double c, const double *nugget_matrix, int n) {
	int count = 0;
	int i = 0;
	while (n > 0) {
		m[i] += c*nugget_matrix[count];
		i+= n; count++;	n--;
	}
}



/*********************************************************************************
logdmvnorm - calculates and returns the log likelihood of observations y given 
	mean matrix mu and var-cov matrix V; -DBL_MAX is returned if the var-cov
	matrix is singular

	V should be in PACKED FORM

	note: matrices y and V will be overwritten
	may want to pass in invV if we have it already
*********************************************************************************/
double logdmvnorm(double *y, const double *mu, double *V, int n) {

	int info = 0;

	double *invV = MATRIX(n,n);
	createIdentityMatrix(invV,n);

	info = LP_sym_pos_solve(V, n, invV, n);  
	if (info != 0) {        // V is singular
		FREE_MATRIX(invV);
		return -DBL_MAX;
	}
	// now invV is inverse of V; V contains cholesky decomp

	double logdet = logDetFromCholesky(V,n);
	double *ans = VECTOR(n);
	vectorSubtract(y, mu, n);  // y is now = y - mu	

	xprimeA(y, invV, ans, n,n); // ans = (y-mu)' %*% invV

	double d = dotprod(ans, y, n);  
	
	FREE_MATRIX(invV);
	FREE_VECTOR(ans);
	return	-(n / 2.0) * (M_LN2 + M_LNPI) - 0.5* (logdet + d);
}



/*********************************************************************************
createCorrMatrix - creates the correlation matrix for design X where

	corr(z(xi), z(xj)) = exp[sum(-B(xi,p - xj,p)**2)], p = 1,2,..,k
	
	The correlation matrix is returned in corr, which must be in packed form 
*********************************************************************************/
void createCorrMatrix(const double *X, const double *B, double *corr, int n, int k) {

	int i, j, p, count = 0;
	zeroPackedMatrix(corr, n);
	double diff;
	double total;
	for (i = 0; i < n; i++) {
	  for(j = i; j < n; j++) {
		diff = total = 0.0;
		for (p = 0; p < k; p++) {
			diff = X[i*k+p] - X[j*k+p];
			total += -B[p]*diff*diff;
		}
		corr[count] = exp(total);
		count++;
	  }
	}
}
#endif


