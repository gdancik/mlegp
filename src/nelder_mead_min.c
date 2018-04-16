/* Nelder-Mead optimization:
   Modified from R's source code (version 2.4.0). Please see at the end of 
   these comments for copyright information.

   Modification Date: 2006/10/21 

   Changes: - removed headers from Rmath.h (and made standalone) 
            - brought in some sort of documentation to help in calling.
            - renamed function to nelder_mead_min to avoid future 
	      potential conflicts in common calls to R.
	    - freed the matrix of vertices P (for some reason, not in the code)
	    
   This program is modified from the source code for the R statistical software
   package, simply so that we can minimize objective functions with external
   fields such as data in them. The modification is minor and listed above.
   	    
   Therefore, the following copyright probably applies. Regardless, the credits
   should go the R team.

 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2006  the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street Fifth Floor, Boston, MA 02110-1301  USA
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "nelder_mead_min.h"
#include "print.h"

#define big 1e+140



#ifdef __useR__
#include <R_ext/Utils.h>
#endif

typedef enum { nmFALSE = 0, nmTRUE } boole;   


void nelder_mead_min(int n, double *Bvec, double *X, double *Fmin,
		     double (*fminfn)(int, double *, void *), 
		     int *fail, double abstol, 
		     double intol, void *ex, double alpha, double bet, 
		     double gamm, int trace, int *fncount, int maxit)
/* int n               = number of parameters
   double *Bvec        = initial values (pointer to 1-D array of p elements)
   double *X           = final value (pointer to 1-D array of p elements)
   double *Fmin        = value at located minimum
   double fminfn       = objective function
   int fail            = 1 if no convergence in maximum number of iterations
   const double abstol = 1e-16;
   const double reltol = 1e-8;
   void *ex            = pointer to external data (perhaps as a struct)
   const double alpha  = 1.0 default  reflection factor
   const double beta   = 0.5 default contraction factor 
   const double gamm   = 2.0 default expansion factor 
   const int trace     = 0; default tracing on if 1.
   int fncount         = number of function evaluations/iterations
   int maxit           = maximum number of iterations  */
{
    char action[50];
    int C;
    boole calcvert;
    double convtol, f;
    int funcount=0, H, i, j, L=0;
    int n1=0;
    double oldsize;
    double **P;
    double size, step, temp, trystep;
    //char tstr[6];
    double VH, VL, VR;

    if (maxit <= 0) {
    	printout("get Fmin..\n");
	*Fmin = fminfn(n, Bvec, ex);
	*fncount = 0;
	*fail = 0;
	printout("\n done, fail = %d\n", fail);
	return;
    }
    if (trace)
     printout("  Nelder-Mead direct search function minimizer\n");
    //MAKE_MATRIX(P,n+1,n+2);
    MAKE_2ARRAY(P, n+1, n+2);

    *fail = nmFALSE;
    f = fminfn(n, Bvec, ex);
    if (!isfinite(f)) {
     printout("oops - function cannot be evaluated at initial parameters\n");
      *fail = nmTRUE;
    } else {
	if (trace)printout("function value for initial parameters = %f\n", f);
	funcount = 1;
	convtol = intol * (fabs(f) + intol);
	if (trace)printout("  Scaled convergence tolerance is %g\n", convtol);
	n1 = n + 1;
	C = n + 2;
	P[n1 - 1][0] = f;
	for (i = 0; i < n; i++)
	    P[i][0] = Bvec[i];

	L = 1;
	size = 0.0;

	step = 0.0;
	for (i = 0; i < n; i++) {
	    if (0.1 * fabs(Bvec[i]) > step)
		step = 0.1 * fabs(Bvec[i]);
	}
	if (step == 0.0) step = 0.1;
	if (trace)printout("Stepsize computed as %f\n", step);
	for (j = 2; j <= n1; j++) {
	    strcpy(action, "BUILD          ");
	    for (i = 0; i < n; i++)
		P[i][j - 1] = Bvec[i];

	    trystep = step;
	    while (P[j - 2][j - 1] == Bvec[j - 2]) {
		P[j - 2][j - 1] = Bvec[j - 2] + trystep;
		trystep *= 10;
	    }
	    size += trystep;
	}
	oldsize = size;
	calcvert = nmTRUE;
	do {
	    #ifdef __useR__
	    R_CheckUserInterrupt();
	    #endif
	    if (calcvert) {
		for (j = 0; j < n1; j++) {
		    if (j + 1 != L) {
			for (i = 0; i < n; i++)
			    Bvec[i] = P[i][j];
			f = fminfn(n, Bvec, ex);
			if (!isfinite(f)) f = big;
			funcount++;
			P[n1 - 1][j] = f;
		    }
		}
		calcvert = nmFALSE;
	    }

	    VL = P[n1 - 1][L - 1];
	    VH = VL;
	    H = L;

	    for (j = 1; j <= n1; j++) {
		if (j != L) {
		    f = P[n1 - 1][j - 1];
		    if (f < VL) {
			L = j;
			VL = f;
		    }
		    if (f > VH) {
			H = j;
			VH = f;
		    }
		}
	    }

	    if (VH <= VL + convtol || VL <= abstol) break;

	    //sprintf(tstr, "%5d", funcount);
	    //if (trace)printout("%s%s %f %f\n", action, tstr, VH, VL);

	    for (i = 0; i < n; i++) {
		temp = -P[i][H - 1];
		for (j = 0; j < n1; j++)
		    temp += P[i][j];
		P[i][C - 1] = temp / n;
	    }
	    for (i = 0; i < n; i++)
		Bvec[i] = (1.0 + alpha) * P[i][C - 1] - alpha * P[i][H - 1];
	    f = fminfn(n, Bvec, ex);
	    if (!isfinite(f)) f = big;
	    funcount++;
	    strcpy(action, "REFLECTION     ");
	    VR = f;
	    if (VR < VL) {
		P[n1 - 1][C - 1] = f;
		for (i = 0; i < n; i++) {
		    f = gamm * Bvec[i] + (1 - gamm) * P[i][C - 1];
		    P[i][C - 1] = Bvec[i];
		    Bvec[i] = f;
		}
		f = fminfn(n, Bvec, ex);
		if (!isfinite(f)) f = big;
		funcount++;
		if (f < VR) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = f;
		    strcpy(action, "EXTENSION      ");
		} else {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = P[i][C - 1];
		    P[n1 - 1][H - 1] = VR;
		}
	    } else {
		strcpy(action, "HI-REDUCTION   ");
		if (VR < VH) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = VR;
		    strcpy(action, "LO-REDUCTION   ");
		}

		for (i = 0; i < n; i++)
		    Bvec[i] = (1 - bet) * P[i][H - 1] + bet * P[i][C - 1];
		f = fminfn(n, Bvec, ex);
		if (!isfinite(f)) f = big;
		funcount++;

		if (f < P[n1 - 1][H - 1]) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = f;
		} else {
		    if (VR >= VH) {
			strcpy(action, "SHRINK         ");
			calcvert = nmTRUE;
			size = 0.0;
			for (j = 0; j < n1; j++) {
			    if (j + 1 != L) {
				for (i = 0; i < n; i++) {
				    P[i][j] = bet * (P[i][j] - P[i][L - 1])
					+ P[i][L - 1];
				    size += fabs(P[i][j] - P[i][L - 1]);
				}
			    }
			}
			if (size < oldsize) {
			    oldsize = size;
			} else {
			    if (trace)
				printout("Polytope size measure not decreased in shrink\n");
			    *fail = 10;
			    break;
			}
		    }
		}
	    }

	} while (funcount <= maxit);

    }

    if (trace) {
      printout("Exiting from Nelder Mead minimizer\n");
      printout("    %d function evaluations used\n", funcount);
    }
    *Fmin = P[n1 - 1][L - 1];
    for (i = 0; i < n; i++) X[i] = P[i][L - 1];
    
	FREE_2ARRAY(P);
    if (funcount > maxit) *fail = 1;
    *fncount = funcount;
}
