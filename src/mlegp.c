/********************************************************
* 12/20/12: all output must be written to R; therefore  *
* I commented all putchars and some of print.h, also 	*
* I changed all exits to returns and commented out all 	*
* asserts						*
********************************************************/

#include <stdio.h>
#include "matrix_vector.h"
#include "gp.h"
#include "fit_gp.h"
#include "print.h"

/******
int main() {

int i;
int n = 5;
double sum;
double *C = PACKED_MATRIX(n);
for (i = 0; i < n*(n+1)/2; i++) {
	C[i] = i;
}

printPackedMatrix("%6.2f ", C, n);

double *reps = VECTOR(n);
for (i = 0; i < n; i++) reps[i] = 1;

//scale_packed_corMatrix(C, 2.0, reps,n);
printout("\n\n");
printVector("%6.2f", reps, n);
printPackedMatrix("%6.2f ", C, n);


return 1;
addNuggetToPackedMatrix(C,10, n);
printf("\n");
printPackedMatrix("%6.2f ", C, n);

double v[] = {1,2,3,4,5};
addNuggetMatrixToPackedMatrix(C,1,v,n);
printPackedMatrix("%6.2f ", C, n);

printf("\n");
int nrowsX = 10, ncolsX = 2;
double vX[] = {1,2,3,4,5,1,2,3,4,5,1,1,1,1,1,2,2,2,2,2};
double *X = MATRIX(nrowsX,ncolsX);
createMatrixByCol(vX, nrowsX,ncolsX,X);
printMatrix("%6.2f ", X, nrowsX,ncolsX, "X=");

int nrowsY = nrowsX;
double Y[] = { 2.3,  5,  6,  7.98,  8,  7.5,  8,  9.023, 10.1, 15.77};

int constantMean = 1, numSimplexTries = 4, maxSimplexIterations = 100;
double simplex_abstol = 1e-16, simplex_reltol = 1e-8; //, *simplex_steps = NULL;
int  BFGS_max_iter = 5;
double BFGS_tol = 0.01, BFGS_h = 0.01; //, BFGS_steps = 0.1;
int rng_seed = 10;

double nugget[1];
nugget[0] = .01;
int nugget_length = 1;
double min_nugget = 0;

int numParams = ncolsX + 1;  // corr params + nugget
int length_estimates = numParams + 2; // add meanReg and sig2
double *estimates = VECTOR(length_estimates);  
int verbose = 2;

int nugget_known = 1;
//int reps = 3;
double MSA = 2.5;
double MSE = 0.5;
fitGP(X, nrowsX, ncolsX, Y, nrowsY,  
	constantMean, 
	numSimplexTries, 
	maxSimplexIterations, 
	simplex_abstol, simplex_reltol,  
	BFGS_max_iter, BFGS_tol, BFGS_h, rng_seed,
	nugget, nugget_length,  min_nugget, estimates, verbose, nugget_known, reps, MSA, MSE);
FREE_MATRIX(C);
FREE_MATRIX(X);
FREE_VECTOR(estimates);

return 0;

}

*****/


