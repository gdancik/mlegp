#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include "print.h"

#define big 1e+140


#define MAKE_2ARRAY(a,m,n) do {                                              \
    size_t ARRAY_H2RESERVED;                                                 \
    MAKE_1ARRAY(a,(m)+1);                                                    \
    if (a==NULL)                                                             \
        break;                                                               \
    (a)[m] = NULL;                                                           \
    for (ARRAY_H2RESERVED=0; ARRAY_H2RESERVED<(m); ARRAY_H2RESERVED++) {     \
        MAKE_1ARRAY((a)[ARRAY_H2RESERVED],(n));                              \
        if ((a)[ARRAY_H2RESERVED]==NULL) {                                   \
            FREE_2ARRAY(a);                                                  \
            break;                                                           \
        }                                                                    \
    }                                                                        \
} while (0)

#define MAKE_1ARRAY(a,n) do {                                                \
    (a) = malloc((n) * sizeof *(a));                                         \
    if ((a)==NULL)                                                           \
        printerr("*** in file %s, function %s(), line %d: "           \
                "out of memory!\n",  __FILE__, __func__, __LINE__);          \
} while (0)                                                                  

#define FREE_1ARRAY(a)  do {                                                 \
    free(a);                                                                 \
    a = NULL;                                                                \
} while (0)


/* ---------- 2D arrays ---------------------- */

#define FREE_2ARRAY(a) do {                                                  \
    size_t ARRAY_H2RESERVED;                                                 \
    if (a==NULL)                                                             \
        break;                                                               \
    for (ARRAY_H2RESERVED=0; (a)[ARRAY_H2RESERVED]!=NULL; ARRAY_H2RESERVED++)\
        FREE_1ARRAY((a)[ARRAY_H2RESERVED]);                                  \
    FREE_1ARRAY(a);                                                          \
} while (0)




void nelder_mead_min(int n, double *Bvec, double *X, double *Fmin,
		     double (*fminfn)(int, double *, void *), 
		     int *fail, double abstol, 
                     double intol, void *ex, double alpha, double bet, 
                     double gamm, int trace, int *fncount, int maxit);

/* int n               = number of parameters
   double *Bvec        = initial values (pointer to 1-D array of n elements)
   double *X           = final value (pointer to 1-D array of n elements)
   double *Fmin        = value at located minimum
   double fminfn       = objective function
   int fail            = 1 if no convergence in maximum number of iterations
   const double abstol = 1e-16;
   const double reltol = 1e-8;
   void *ex            = pointer to external data (perhaps as a struct)
   const double alpha  = 1.0 default,  reflection factor
   const double beta   = 0.5 default, contraction factor 
   const double gamm   = 2.0 default, expansion factor 
   const int trace     = 0; default, tracing on if 1.
   int fncount         = number of function evaluations/iterations
   int maxit           = maximum number of iterations  */
