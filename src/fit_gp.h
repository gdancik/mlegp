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

/***************************************************************************
fit_gp.h - fits a gp to a set of observations, using numerical methods and 
	closed form solutions to obtain mles
***************************************************************************/

#ifndef __FIT_GP__
#define __FIT_GP__

#include <time.h>
#include <math.h>

#ifdef __useR__
#include <R_ext/Utils.h>
#endif

#include "matrix_vector.h"
#include "gp.h"
#include "print.h"
#include "nelder_mead_min.h"
#include "SFMT.h"
#include "lbfgs.h"

// parameters that do not vary; will be passed into minimization functions (f_min, etc)
typedef struct {
	double *X;       	// the design that produces the observations 
	double *fX;      	// for the mean function will have the form (fX)'B
	double *Y;		// the observations
	int nY;                 // number of observations (should be mean observations in the case of reps)
	int ncolsX;		// number of columns of X 
	int ncolsfX; 		// number of columns of fX
	double h;		// used in derivative approximation [f(x+h) - f(x+h)] / h
	double *nuggetMatrix;   // pointer to nugget 'matrix'; NULL otherwise
	int verbose;		// level of detail for output
	// estimateNugget = 1 if we are estimating a nugget term (either constant or scalar multiplier of nug matrix)
	// to specify a known nugget, assign a nuggetMatrix but set estimateNugget to 0	
	int estimateNugget;   
	int nugget_known;
	double min_nugget;	// minimum value of the nugget, or can be used to force a nugget
	// if we know the nugget!
	// # of reps per design point (treat as double to avoid silly calculation mistakes)
	//double *reps;		
	//double MSA;
	//double MSE;
} gp_params; 


void gp_params_summary(gp_params *p) {
	printout("summary of gp_params object\n");
	printout("num obs: %d, ncolX: %d, ncolfX: %d\n", p->nY, p->ncolsX, p->ncolsfX);
	printMatrix("%6.3f ", p->X, p->nY, p->ncolsX, "X=");
	printMatrix("%6.3f ", p->Y, p->nY, 1, "Y=");
	printMatrix("%6.3f ", p->fX, p->nY, p->ncolsfX, "fX=");
}


/* vector_exp_check - takes exp of the first n elements in v; sets to 0 if element < -500 */
void vector_exp_check(double *v, int n) {
	int i = 0;
	for (i = 0; i < n; i++) {
		if (v[i] < -500) v[i] = 0;
		else v[i] = exp(v[i]);
	}
}

/* vector_log - takes log of the first n elements in v */
void vector_log(double *v, int n) {
	int i = 0; 
	for (i = 0; i < n; i++) {
		v[i] = log(v[i]);
	}
}


/**********************************************************************************
f_min - finds the neg log like of the gp
	const int num_p = the # of parameters to estimate where closed form mle's
			  do not exist (always equal to ncolX + 1(for a nugget term))
	double *orig_v - vector of parameters being estimated, consisting of
		correlation coefficients, sig2(nugget), in that order; all 
		parameters must be on log scale; sig2(nugget) is a scaled nugget 
		and is actually sig2(nugget) / sig2(gp) 
	void *params - pointer to gp_params structure, which determines whether
			we are using a constant nugget, nugget matrix, the
			mean function, etc.
	
	returns: -log likelihood of gp; or DBL_MAX if there is an error inverting 
			var-cov matrix
***********************************************************************************/
double f_min(const int num_p, double *orig_v, void *params) {
	int i = 0;
	gp_params *p = (gp_params *) params;

//	gp_params_summary(p);

	double *corr = PACKED_MATRIX(p->nY);

	double *v = VECTOR(num_p);

	vectorCopy(orig_v, v, num_p);
	vector_exp_check(v, num_p);   // transform back to original values

	double *B = VECTOR(p->ncolsX);
	for (i = 0; i < p->ncolsX; i++) {
		B[i] = v[i];
	}
	
	createCorrMatrix(p->X, B, corr, p->nY, p->ncolsX);

	if (p->estimateNugget == 1) {  // nugget is specified up to a constant
		if (p->nuggetMatrix != NULL) addNuggetMatrixToPackedMatrix(corr, v[p->ncolsX], p->nuggetMatrix, p->nY);
		else if (p->estimateNugget == 1) addNuggetToPackedMatrix(corr, v[p->ncolsX], p->nY);
	}
	else if (p->estimateNugget == 0 && p->nuggetMatrix != NULL) {  // nugget is KNOWN exactly
		double sig2 = 1;
		if (p->nugget_known == 9) {
			printerr("calcANOVAsig2 not implemented\n");
			return(-1); //exit(-1);
			//sig2 = calcANOVASig2(p->MSA, p->MSE, corr, p->nY, p->reps[1]); // make sure corr is really correlations only
		}
		else if (p->nugget_known == 1) sig2 = v[p->ncolsX];
		if (!isfinite(sig2)) {   
			FREE_VECTOR(B);
			FREE_VECTOR(v);
			FREE_MATRIX(corr);
			return DBL_MAX;
		}
		packed_matrix_scale_const(corr, sig2, p->nY);
		addNuggetMatrixToPackedMatrix(corr, 1.0, p->nuggetMatrix, p->nY);
	}
		addNuggetToPackedMatrix(corr, p->min_nugget, p->nY); // add the minimum nugget

	double *corrInv = MATRIX(p->nY, p->nY);
	createIdentityMatrix(corrInv, p->nY);

	double *copyCorr = PACKED_MATRIX(p->nY);
	copyPackedMatrix(corr, copyCorr, p->nY);

	int info = LP_sym_pos_solve(copyCorr, p->nY, corrInv, p->nY);
	if (info != 0) {  // copyCorr is singular 
		FREE_VECTOR(B);
		FREE_VECTOR(v);
		FREE_MATRIX(corr);
		FREE_MATRIX(copyCorr);
		FREE_MATRIX(corrInv);
		return DBL_MAX;
	}
	// copyCorr contains cholesky decomp, corrInv is corrInv (or inverse of covariance matrix) 
	
	double *bhat = VECTOR(p->ncolsfX);	

	if (calcBhat(p->fX, p->nY, p->ncolsfX, corrInv, p->Y, bhat) != 0) { 
		FREE_VECTOR(B);
		FREE_VECTOR(v);
		FREE_MATRIX(corr);
		FREE_MATRIX(copyCorr);
		FREE_MATRIX(corrInv);
		return DBL_MAX;
	}
	double *mu_hat = MATRIX(p->nY, 1);
	matrix_multiply(p->fX, p->nY, p->ncolsfX, bhat, 1, mu_hat);

	if (p->nugget_known == 0) {	
		double sig2 = calcMLESig2(p->Y, mu_hat, corrInv, p->nY);
		packed_matrix_scale_const(corr, sig2, p->nY);
	}
	double *Y = VECTOR(p->nY);
	copyVector(p->Y, Y, p->nY);

	//printout("get ans = ");
	double ans = - logdmvnorm(Y, mu_hat, corr, p->nY);
	//printout("%f\n", ans);
	FREE_MATRIX(corr);
	FREE_VECTOR(B);
	FREE_VECTOR(v);
	FREE_MATRIX(corrInv);
	FREE_MATRIX(copyCorr);
	FREE_VECTOR(bhat);
	FREE_MATRIX(mu_hat);
	FREE_VECTOR(Y);

	return ans;

} // end f_min


// for printing progress after each iteration in the BFGS method; returning anything other than a 0
// will cancel the optimization process
static int BFGS_progress(
	void *instance,		       // pointer to external data
	const lbfgsfloatval_t *x,      // parameter vector
	const lbfgsfloatval_t *g,      // the gradient vector
	const lbfgsfloatval_t fx,      // function value
	const lbfgsfloatval_t xnorm,   // euclidian norm of variables
	const lbfgsfloatval_t gnorm,   // euclidian norm of gradient
	const lbfgsfloatval_t step,    // line search step size
	int n,                         // # of parameters
	int k,                         // iteration count
	int ls                         // # of evaluations called for this iteration
	)
{
	printout("\titeration: %d,", k);
	printout("loglike = %f\n", -fx);
	return 0;
}


/***********************************************************************************
fdf_evaluate - used with BFGS method; calculates gradient vector and evaluates
 	       the function at the current parameter value; partial derivatives
	       are approximated using a forward difference approximation for each
               parameter. If forward difference approximation cannot be evaluated,
               a backward difference approximation is tried; if none of these work
               or if the current parameter vector cannot be evaluated, the gradient
               is set to zero. This function is not called until after simplex method, 
	       and setting the gradient to 0 is consistent with having an mle near the edge of 
	       a parameter space where the var-cov matrix becomes singular.
************************************************************************************/
static lbfgsfloatval_t fdf_evaluate(
	void *params,              // pointer to gp_param struct
	const lbfgsfloatval_t *v,  // current values (on log scale)
	lbfgsfloatval_t *g,        // gradient vector
	const int n,               // number of parameters
	const lbfgsfloatval_t step // current step of line search (not used by me)
	)
{

	double fv = f_min(n, (double* )v, params);          // current value

	double *v_plus_h = VECTOR(n);
	gp_params *p = (gp_params *) params;

	int i;
	double fv_plus_h; 
	double *v_copy = VECTOR(n);
	copyVector(v, v_copy,n);         // on log scale
	
	vector_exp_check(v_copy, n);      // v_copy is on actual scale
	for (i = 0; i < n; i++) {
		copyVector(v_copy, v_plus_h, n);   // v_plus_h contains values on actual scale
		v_plus_h[i] = v_copy[i] + p->h;
		vector_log(v_plus_h, n);           // put back on log scale
		fv_plus_h = f_min(n, v_plus_h, params);
		if (fv == DBL_MAX) {
			g[i] = 0;
		}
		else if (fv_plus_h == DBL_MAX) {     // try backward approximation
			copyVector(v_copy, v_plus_h, n);   // v_plus_h contains values on actual scale
			v_plus_h[i] = v_copy[i] - p->h;
			vector_log(v_plus_h, n); // put back on log scale
			fv_plus_h = f_min(n, v_plus_h, params);
			if (fv_plus_h == DBL_MAX) {
			    g[i] = 0;
			}
			else g[i] = (fv - fv_plus_h) / -p->h;		

		}
		else g[i] = (fv_plus_h - fv) / p->h;
	}

	FREE_VECTOR(v_plus_h);
	FREE_VECTOR(v_copy);
	return fv;
}


/* findMinEuclidianDist - finds the minimum (and max) Euclidian distance between all pairs of points
	in the design matrix X (excluding identical design points). After returning, m1 points 
        to the minimum distance and m2 points to the maximum distance
*/
void findMinEuclidianDist(const double *X, const int nrows, const int ncols, double *m1, double *m2) {
	int i, j, p;
        *m1 = DBL_MAX;
        *m2 = 0;
	double d, diff;
        for (i = 0; i < nrows - 1; i++) {
                for (j = i+1; j < nrows; j++) {
			d = 0.0;
			for (p = 0; p < ncols; p++) {
				diff = MATRIX_GET(X,i,p,ncols) - MATRIX_GET(X,j,p,ncols);
				d+= diff*diff;
			}
                        if (d > 0 && d < *m1) {
                                *m1 = d;
                        }
                        if (d > *m2) {
                                *m2 = d;
                        }
                }
        }
}

/****************************************************************************************************
getUnivariateCorRange
	const double *X - pointer to design matrix 
	double *m1 - at end of function call, pointer to value of beta for which max correlation btwn
			any 2 observations is 0.65
	double *m2 - at end of function call, pointer to value of beta for which min correlation btwn
			any 2 observations is 0.3

	note: the purpose of this function is to find initial values of correlation parameters so 
		that the var-cov matrix is not (approximately) singular (which could happen if correlation btwn any
		2 observations is close to 1)
****************************************************************************************************/
void getUnivariateCorRange(const double *X, int nrows, int ncols, double *m1, double *m2) {
   double xmin, xmax;
   findMinEuclidianDist(X, nrows, ncols, &xmin, &xmax);
   *m1 = -log(.65) / xmin;    
   *m2 = -log(.3) / xmin;     
}


void printBFGSReturnMsg(int x) {

switch(x) {

	case 0:
	case LBFGSERR_ROUNDING_ERROR:           // most likely this is convergence in line search
	case LBFGSERR_MINIMUMSTEP:              // line search step became smaller than lbfgs_parameter_t::min_step
        case LBFGSERR_MAXIMUMSTEP:               // line search step became larger than lbfgs_parameter_t::max_step
		return;
        default: 
		printout("L-BFGS warning: ");
		if (x == LBFGSERR_LOGICERROR) printerr("logic error\n");
		else if (x == LBFGSERR_OUTOFMEMORY) printerr("ERROR: out of memory\n");
		else if (x == LBFGSERR_MAXIMUMLINESEARCH) printout("line search routine reached max # of evaluations\n");
		else if (x == LBFGSERR_WIDTHTOOSMALL) 
			printout("relative width of interval of uncertainty is at most lbfgs_parameter_t::xtol\n");
		/************** check input parameters***************************************************************/
		else if (x == LBFGSERR_INVALID_N) printerr("invalid number of variables specified\n");
		else if (x == LBFGSERR_INVALID_N_SSE) printerr("invalid number of variables (for SSE) specified\n"); 
		else if (x == LBFGSERR_INVALID_MINSTEP) 
			printerr("Invalid parameter lbfgs_parameter_t:min_step specified\n");
		else if (x == LBFGSERR_INVALID_MAXSTEP) 
			printerr("invalid parameter lbfgs_parmaeter_t:max_step specified\n");
		else if (x == LBFGSERR_INVALID_FTOL) 
			printerr("invalid parameter lbfgs_parameter_t::ftol specified\n");
		else if (x == LBFGSERR_INVALID_GTOL) 
			printerr("invalid parameter lbfgs_parameter_t::gtol specified\n");
		else if (x == LBFGSERR_INVALID_XTOL) 
			printerr("invalid parameter lbfgs_parameter_t::xtol specified\n");
		else if (x == LBFGSERR_INVALID_MAXLINESEARCH) printerr("lbfgs_parameter_t::max_linesearch specified\n");
		else if (x == LBFGSERR_INVALID_ORTHANTWISE) 
			printerr("invalid parameter lbfgs_parameter_t::orthantwise_c specified\n");
		else if (x == LBFGSERR_INVALIDPARAMETERS) printerr("logic error (negative line-search step) occurred\n");
		/****************************************************************************************************/
		else if (x == LBFGSERR_OUTOFINTERVAL) printout("line search step went out of interval of uncertainty\n");
		else if (x == LBFGSERR_INCORRECT_TMINMAX) 
			printout("logic error occured or interval of uncertainty became too small\n");
		else if (x == LBFGSERR_MAXIMUMITERATION) 
			printout("algorithm reaches maximum # of iterations\n");
		else if (x == LBFGSERR_INCREASEGRADIENT) 
			printout("current search direction increases object function value\n");
		else printerr("warning: unknown error message returned from lbfgs\n");
	}
	/************************ NOT USED ********************
	else if (x == LBFGSERR_CANCELED) printerr("canceled by user\n");
	else if (x == LBFGSFALSE) printout("FALSEi\n");
	else if (x == LBFGSTRUE) printout("TRUE\n");
	else if (x == LBFGSERR_UNKNOWNERROR) printerr("unknown error\n");
	*****************************/
}


/***************************************************************************************************************
fitGP - the driver function that finds mles for gp parameters; we run 'numSimplexTries' simplexes to find a good initial value, then continue with bfgs method
	const double *X - pointer to design matrix for the respone(s)
	const int nrowsX - # of rows of X, = to number of observations in Y 
	const int ncolsX - # of columns of X, = to number of design parameters	
	const double *Y - pointer to 1-column matrix of observations
	const int nrowsY - number of observations (must match nrowsX)
	int constantMean - 1 if gp should have a constant mean; otherwise mean function will be (1 X) %*% B
		(linear regression in X with intercept added)
	int numSimplexTries - the number of simplexes to run before moving to bfgs method
	int maxSimplexIterations - maximum number of iterations / simplex
	double simplex_abstol - absolute tolerance for simplex (will be set to -DBL_MAX)
	double simplex_reltol - relative tolerance of simplex
	int BFGS_max_iter - max # of iterations for bfgs method
	double BFGS_tol - stopping criteria for bfgs method; minimization terminates when 
	 			||g|| <  BFGS_tol * max(1, ||x||), where ||.|| denotes the 
   				Euclidean (L2) norm. 
	double BFGS_h -   for finite difference approximation to calculate gradient during BFGS
	int rng_seed - random # seed
	double *nugget - pointer to initial nugget value or nugget matrix; this MUST be a nugget matrix if nugget_known is 1 
	int nugget_length - 1 for constant nugget, or greater than 1 for nugget matrix, in which case
		 should be equal to the # of observations; 0 otherwise
	double min_nugget - minimum value of nugget (added to diagonal of var-cov matrix)
	double *estimates - stores values of mles at end of function call; order of estimates will be:
			    ( meanReg parameters, correlation parameters, sig2(GP), sig2(nugget) )
	int verbose - 0 for no status reporting; 1 for some; 2 for most
	int nugget_known = 1 if nugget is fixed and to find MLE of sig2; or 9 if nugget is fixed and to find MOM estimate of Sig2 (no longer implemented)
	returns 0 if success; -1 otherwise

**************************************************************************************************************/
int fitGP(const double *X, const int nrowsX, const int ncolsX, 
	   const double *Y, int nrowsY, 
	   int constantMean, 
	   int numSimplexTries, int maxSimplexIterations, double simplex_abstol, double simplex_reltol, 
	   int BFGS_max_iter, double BFGS_tol, double BFGS_h,  
           int rng_seed, 
	   double *nugget, int nugget_length, double min_nugget, 
	   double *estimates, int verbose, int nugget_known/*, double *reps, double MSA, double MSE*/) {

	init_gen_rand(rng_seed);
	simplex_abstol = -DBL_MAX;

	if (nrowsX != nrowsY) {
		printerr("ERROR: number of observations does not match number of design points\n");
		return -1;
	}
	if (nugget_length > 1 && nugget_length != nrowsY) {
		printerr("ERROR: length of nugget matrix must match number of observations, or use constant nugget\n");
		return -1;
	}

	double *fX;
	int ncolsfX = 0;
	if (constantMean > 0) {
		fX = MATRIX(nrowsX, 1);
		ncolsfX = 1;
		setMatrix(fX, nrowsX, ncolsfX, 1.0);
	}
	else {
		fX = MATRIX(nrowsX, ncolsX+1);
		ncolsfX = ncolsX+1;
		double *ones = MATRIX(nrowsX, 1);
		setMatrix(ones, nrowsX, 1, 1.0);
		cbind(ones, X, fX, 1, ncolsX, nrowsX);
		FREE_MATRIX(ones);
	}

	gp_params gp;
	gp.X = (double *) X;	   // design matrix for observed responses
	gp.nY = nrowsX;
	gp.ncolsX = ncolsX;
	gp.ncolsfX = ncolsfX;
	gp.fX = (double *) fX;	           // design matrix for mean
	gp.Y = (double *) Y;
	gp.ncolsX = ncolsX;  // number of correlation parameters 
	gp.h = BFGS_h;       // for finite difference approximation
	gp.nuggetMatrix = NULL;    // pointer to nugget array for non-constant nugget 
	if (nugget_length > 1) {
		gp.nuggetMatrix = nugget;
		//return -1;
	}
	gp.verbose = verbose;
	gp.estimateNugget = 1;
	if (nugget_length < 1) gp.estimateNugget = 0;  // we are NOT estimating nugget
	if ( (nugget_known == 1) || (nugget_known == 9)) gp.estimateNugget = 0;  // we are NOT estimating nugget
	gp.nugget_known = nugget_known;
	gp.min_nugget = min_nugget;

	//gp.reps = reps;
	//gp.MSA = MSA;
	//gp.MSE = MSE;

	//printout("nugget: ");
	//printVector("%2.2f ", gp.nuggetMatrix, gp.nY);
	//printout("reps: ");
	//printVector("%2.1d ", (double*) gp.reps, gp.nY);
	//for (int i = 0; i < gp.nY; i++)
	//	printout("%d ", gp.reps[i]);
	//printout("\n");	
	//exit(-1);

	// the vector of parameters to maximize; correlation parameters + nugget 
	double *v;
	int v_length = ncolsX;       // for correlation parameter estimates

	if ((gp.estimateNugget == 1) || (gp.nugget_known == 1)) v_length++; // either for nugget or sig2   

	v = VECTOR(v_length);  

	double initial_scaled_nugget = 0;       // we will use the same initial nugget for each simplex
	if (gp.estimateNugget == 1) {
		// choose initial values 
		double sig2 =  vectorVariance(Y, nrowsY);   
		initial_scaled_nugget = 0;       // we will use the same initial nugget for each simplex

		if (nugget_length == 1 && nugget[0] > 0) {   // we are adding a constant nugget
			double nug;
			if (nugget_length == 1) nug =  *nugget; 
			else nug = 1.0;
			nug = nug / sig2;	// scale the nugget
			v[ncolsX] = nug;         
			initial_scaled_nugget = v[ncolsX];
		}
		else if (nugget_length > 1) {
			v[ncolsX] = 1.0 / sig2;    // scaled nugget is sig2e / sig2GP
			initial_scaled_nugget = v[ncolsX];
		}	
		printout("intial_scaled nugget is %f\n", v[ncolsX]);
	}

	double m1 = 0, m2;
	// correlation btwn 2 closest points is btwn m1 and m2, assuming a constant beta 
	getUnivariateCorRange(X, nrowsX, ncolsX, &m1, &m2);   

	int simplex_try;
	double fval = 0;
	double best_fval = DBL_MAX;
	int best_try = 0;
	double *best_v = VECTOR(v_length);
	double *best_v_from_simplex = VECTOR(v_length);
	double *v_save = VECTOR(v_length);
	int i;

  	const double alpha  = 1.0; // default  reflection factor for nelder-mead 
	const double bet   = 0.5; // default contraction factor for nelder-mead
	const double gamm  = 2.0; // default expansion factor for nelder-mead
	int maxit = maxSimplexIterations;

if (1) {   /// MAXIMIZATION PROCEDURE	
	for (simplex_try = 0; simplex_try < numSimplexTries; simplex_try++) { 
	
		if (verbose > 0 && maxit > 0) {
			printout("running simplex # %d...\n", simplex_try+1);
		}

		//printout("select initial Beta between %f, %f\n", m1,m2);
		for (i = 0; i < ncolsX; i++) {
			v[i] = m1 + (m2 - m1) * genrand_res53(); 
		}
	
		if (v_length > ncolsX) {     // we are estimating a nugget term or sig2

			if(gp.nugget_known == 0) v[ncolsX] = initial_scaled_nugget;   
			if(gp.nugget_known == 1) { // initial sig2
				v[ncolsX] = vectorVariance(Y, nrowsY);   
			}
		}
		
	//	printout("The %d initial params are\n", v_length);
	//	printVector("%2.3f ", v, v_length);

		// put all initial param values on the log scale
		vector_log(v, v_length);

		//printout("On log scale:\n");
		//printVector("%2.3f", v, v_length);
	
		//printout("evaluate...\n");
		//double ans = -f_min(v_length, v, &gp);
		//printout("...done, -loglike = %f\n", ans);

		int fail = 0;
		
		int trace = 0;   // trace off
		if (verbose == 2) trace = 1; 

		int fncount;

	        if (maxit > 0) {
			//printout("call nelder_mead_min...\n");
		  nelder_mead_min(v_length, v, v_save, &fval,
		     f_min, &fail, simplex_abstol, simplex_reltol, &gp, alpha, bet, 
		     gamm, trace, &fncount, maxit);
			printout("...done\n");
			
		  if (verbose > 0) {
		    	if (fail != 1) 
				printout("...simplex #%d complete, loglike = %f (convergence)\n", simplex_try+1, -fval);
		        else printout("...simplex #%d complete, loglike = %f\n", simplex_try+1, -fval);
		  }	
		}

		if (maxit == 0) vectorCopy(v, v_save, v_length);   // because v_save gets copied to best_v	
		if (simplex_try == 0 || fval < best_fval) {
			best_fval = fval;
			best_try = simplex_try + 1;
			vectorCopy(v_save, best_v, v_length);
		}
		if (maxit == 0) break;       // initial values must be in best_fval		
	}

	vectorCopy(best_v, best_v_from_simplex, v_length);
	if (verbose > 0 && maxit > 0 ) {
		printout("\nusing L-BFGS method from simplex #%d...\n", best_try);
	}

	lbfgs_parameter_t *BFGS_params = malloc(sizeof(lbfgs_parameter_t));
	lbfgs_parameter_init(BFGS_params);

	BFGS_params->epsilon = BFGS_tol; //1e-5;
	BFGS_params->max_iterations = BFGS_max_iter; //note 0, leads to iteration until convergence (or error)

	int ret;
	if (BFGS_params->max_iterations > 0) {
		if (verbose == 0) ret = lbfgs(v_length, best_v, fdf_evaluate, NULL, (void *) &gp, BFGS_params);
		else {
			ret = lbfgs(v_length, best_v, fdf_evaluate, BFGS_progress, (void *) &gp, BFGS_params);
			/* Report the result. */
			printout("...L-BFGS method complete\n");
			printBFGSReturnMsg(ret);
		}
	}

	///// REPORT RESULTS ////
	fval = f_min(v_length, best_v, (void *) &gp);  // do this while estimates are still on log scale
	if (isnan(fval) || fval == DBL_MAX) {      // if lbfgs method does is not successful, use nelder mead result
		vectorCopy(best_v_from_simplex, best_v, v_length);
		fval = f_min(v_length, best_v, (void *) &gp);  // do this while estimates are still on log scale
	        printout("warning: log likelihood in L-BFGS method returns NaN; estimates from simplex method will be used\n"); 
	}

	if (verbose) printout("\nMaximum likelihood estimates found, log like =  %f\n", -fval);   

	// put best estimates on actual scale
	vector_exp_check(best_v, v_length);
}/// end maximization

	// get mle's for beta corr parameters	
	double *B = VECTOR(ncolsX);

	for (i = 0; i < ncolsX; i++) {
		B[i] = best_v[i];
	}	

//	printout("now B = ");
//	printVector("%f", B, ncolsX);
	double *corr = PACKED_MATRIX(nrowsX);	
	createCorrMatrix(X, B, corr, nrowsX, ncolsX); 

	double sig2 = -1.0;
	if (gp.estimateNugget == 1) {
		printout("addNuggets...\n");
   	  	if (gp.nuggetMatrix != NULL) addNuggetMatrixToPackedMatrix(corr, best_v[gp.ncolsX], gp.nuggetMatrix, gp.nY);
          	else if (gp.estimateNugget == 1) addNuggetToPackedMatrix(corr, best_v[gp.ncolsX], gp.nY);
	}
	else if (gp.estimateNugget == 0 && gp.nuggetMatrix != NULL) {  // nugget is KNOWN
		sig2 = -1.0;
		if (gp.nugget_known == 9) {
			printerr("calcANOVASig2 not implemented\n");
			return(-1); //exit(-1);
			//sig2 = calcANOVASig2(gp.MSA, gp.MSE, corr, gp.nY, gp.reps[1]); // make sure corr is really correlations only
		}
		else if (gp.nugget_known == 1) {
			sig2 = best_v[gp.ncolsX];
		}

		if (!isfinite(sig2)) {   
			printerr("ERROR: final sig2 is not finite!\n");
			FREE_VECTOR(B);
			FREE_VECTOR(v);
			FREE_VECTOR(best_v);
			FREE_VECTOR(best_v_from_simplex);
			FREE_VECTOR(v_save);
			FREE_MATRIX(fX);
			FREE_MATRIX(corr);
			return -1;
		}

		packed_matrix_scale_const(corr, sig2, gp.nY);
		addNuggetMatrixToPackedMatrix(corr, 1.0, gp.nuggetMatrix, gp.nY);
	}
        addNuggetToPackedMatrix(corr, gp.min_nugget, gp.nY); // add the minimum nugget

	double *corrInv = MATRIX(nrowsX, nrowsX);
	createIdentityMatrix(corrInv, nrowsX);
	double *copyCorr = PACKED_MATRIX(nrowsX);
	copyPackedMatrix(corr, copyCorr, nrowsX);

	int info = LP_sym_pos_solve(copyCorr, nrowsX, corrInv, nrowsX);
	if (info != 0) {  // singular or non-pos matrix!
		printerr("ERROR: final var-cov matrix is singular!\n");
		FREE_VECTOR(B);
		FREE_VECTOR(v);
		FREE_VECTOR(best_v);
		FREE_VECTOR(best_v_from_simplex);
		FREE_VECTOR(v_save);
		FREE_MATRIX(fX);
		FREE_MATRIX(corr);
		FREE_MATRIX(corrInv);
		FREE_MATRIX(copyCorr);
		return -1;
	}
	// copyCorr contains cholesky decomp, corrInv is corrInv
	
	double *bhat = VECTOR(ncolsfX);

	if (calcBhat(fX, nrowsX, ncolsfX, corrInv, Y, bhat) != 0) {
		FREE_VECTOR(B);
		FREE_VECTOR(v);
		FREE_VECTOR(best_v);
		FREE_VECTOR(best_v_from_simplex);
		FREE_VECTOR(v_save);
		FREE_MATRIX(fX);
		FREE_MATRIX(corr);
		FREE_MATRIX(corrInv);
		FREE_MATRIX(copyCorr);
		FREE_VECTOR(bhat);
		printerr("ERROR: final scaled var-cov matrix is singular in calcBhat!\n");
		return -1;
	}	

	double *mu_hat = MATRIX(nrowsX, 1);
	matrix_multiply(fX, nrowsX, ncolsfX, bhat, 1, mu_hat);

	if (gp.nugget_known == 0) {
		sig2 = calcMLESig2(Y, mu_hat, corrInv, nrowsX);	
	}
	// create solution vector
	int index = 0;
	for (i = 0; i < ncolsfX; i++) {
		estimates[index] = bhat[i];   // meanReg params
		index++;
	}
	for (i = 0; i < ncolsX; i++) {
		estimates[index] = B[i];     // corr params
		index++;
	}
	estimates[index] = sig2;                   // sig2(GP)
	index++;
	// if necessary, add the true (unscaled) value of the nugget to the estimates list
	if (gp.estimateNugget == 1) {
		estimates[index] = (best_v[ncolsX]) * sig2;
		index++;
	}

	FREE_VECTOR(B);
	FREE_VECTOR(v);
	FREE_VECTOR(best_v);
	FREE_VECTOR(best_v_from_simplex);
	FREE_VECTOR(v_save); 
	FREE_VECTOR(bhat);
	FREE_MATRIX(fX);
	FREE_MATRIX(corr);
	FREE_MATRIX(corrInv);
	FREE_MATRIX(copyCorr);
	FREE_MATRIX(mu_hat);
	return 0;

}

void fitGPTest(double *in_X, double *in_Y) {
	printout("fitGPTest, now exit..\n");
}


/*
fitGPfromR - the function actually called from R, that calls fitGP. See fitGP for parameter
		descriptions
*/
void fitGPfromR(double *X, int *nrowsX, int *ncolsX, 
	   double *Y, int *nrowsY, 
	   int *constantMean, 
	   int *numSimplexTries, int *maxSimplexIterations, double *simplex_abstol, double *simplex_reltol, 
	   int *BFGS_max_iter, double *BFGS_tol, double *BFGS_h,  
           int *rng_seed, 
	   double *nugget, int *nugget_length, double *min_nugget, 
	   double *estimates, int *verbose, int *nugget_known, /*double *reps, double *MSA, double *MSE,*/ int *ret) {

	   double *X_mat = MATRIX(*nrowsX, *ncolsX);
	   createMatrixByCol(X, *nrowsX, *ncolsX, X_mat);

	   *ret = fitGP(X_mat, *nrowsX, *ncolsX, Y, *nrowsY, *constantMean, 
		 *numSimplexTries, *maxSimplexIterations, *simplex_abstol, *simplex_reltol,
		 *BFGS_max_iter, *BFGS_tol, *BFGS_h,
		 *rng_seed,
		 nugget, *nugget_length, *min_nugget, estimates, *verbose, *nugget_known/*, reps, *MSA, *MSE*/);

	   FREE_MATRIX(X_mat); 
}

#endif


