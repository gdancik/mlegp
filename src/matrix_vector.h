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
matrix.h - provides several matrix and vector functions, some of which use 
	lapack and cblas
***************************************************************************/

#ifndef __MATRIX_VECTOR__
#define __MATRIX_VECTOR__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "array.h"
#include "lapack_headers.h"
#include "print.h"

#define MATRIX(a,b) (double *) malloc((a)*(b)*sizeof(double))

// stores symmetric matrix in packed format; must be used for inversion / det of var-cov matrix
#define PACKED_MATRIX(n) (double *) malloc((n) * ((n) +1)/2 *sizeof(double))
#define VECTOR(n) (double *) malloc((n)*sizeof(double))

#define MATRIX_GET(m,i,j,ncol) m[(i)*(ncol)+(j)]
#define MATRIX_SET(m,i,j,ncol,a) m[(i)*(ncol)+(j)] = (a)

#define FREE_MATRIX(m) free(m) 
#define FREE_VECTOR(v) free(v) 



/************************************************************************
createMatrixByRow -  
	d - array of elements in ROW order
	nrow - # of rows
	ncol - # of columns
	m - pointer to a matrix that will store the new matrix, 
************************************************************************/
void createMatrixByRow(const double *d, int nrow, int ncol, double *m) {
	int i, j;
	int count = 0;
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			m[i*ncol + j] = d[count];
			count++;
		}
	}
}
/************************************************************************
createMatrixByCol -  
	d - array of elements in COLUMN order
	nrow - # of rows
	ncol - # of columns
	m - pointer to a matrix that will store the final matrix, 
************************************************************************/
void createMatrixByCol(const double *d, int nrow, int ncol, double *m) {
	int i, j;
	int count = 0;
	for (i = 0; i < ncol; i++) {
		for (j = 0; j < nrow; j++) {
			m[j*ncol + i] = d[count];
			count++;
		}
	}
}


/************************************************************************
createPackedMatrix -  
	d - vector of elements to store in packed matrix format
	n - # of rows/columns
	m - pointer to a matrix that will store the final matrix 
************************************************************************/
void createPackedMatrix(const double *d, int n, double *m) {
	int i;
	for (i = 0; i < n*(n+1)/2; i++) {
			m[i] = d[i];
	}
}

void zeroPackedMatrix(double *m, int n) {
	int i;
	for (i = 0; i < n*(n+1)/2; i++) {
		m[i] = 0;
	}
}

double sumPackedMatrix(const double *m, const int n) {
	int i = 0, count = 0;
	double sum = 0.0, diag = 0.0;
	for (i = 0; i < n*(n+1)/2; i++) {
		sum += m[i];
	}
	
	// get diag
	for (i = n; i > 0; i--) {
		diag += m[count];
		count+=i;
	}

	return (sum - diag)*2.0 + diag;	
}

void createIdentityMatrix(double *m, int n) {
	int i,j;
	for (i = 0; i < n; i++) {
		for(j = 0; j< n; j++){
			if (i == j) m[i*n+j] = 1;
			else m[i*n+j] = 0;
		}
	}
}

// sets all elements of matrix m to value
void setMatrix(double *m, const int n1, const int n2, double value) {
	int i;
	for (i=0; i < n1*n2;i++) {
		m[i] = value;
	}
}

// pack and unpack only for symmetric matrices
void unpackMatrix(const double *a, double *b, int n) {
	int i, j, count=0;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			b[i*n+j] = b[j*n+i] = a[count];
			count++;	
		}
	}
}
// will contain upper right of matrix, from left to right
void packMatrix(const double *a, double *b, int n) {
	int i, j, count = 0; 
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			b[count] =  a[i*n+j];
			count++;
		}
	}
}

// copy packed matrix m1 to packed matrix m2
void copyPackedMatrix(const double *m1, double *m2, int n) {
		int i;
		for (i = 0; i < n*(n+1)/2; i++) {
			m2[i] = m1[i];
		}
}

void copyVector(const double *v1, double *v2, int n) {
	int i;
	for (i = 0; i < n; i++) {
		v2[i] = v1[i];
	}
}

/*************************************************************************
LP_sym_pos_solve - solves AX=B for a packed (symmetric) matrix A
	A - pointer to a PACKED matrix
	B - a (unpacked) matrix to store the final answer
	n - number of rows/columns of a  
	nrhs - number of columns of b

	returns 1 if A cannot be inverted; (probably) 0 otherwise
*************************************************************************/
int LP_sym_pos_solve(double *A, int n, double *B, int nrhs)
{
	char uplo = 'L'; // L
	int info;
	dppsv_(&uplo, &n, &nrhs, A, B, &n, &info);
	return info;
}


/* C front end to LAPACK's DGESV routine.
 * Computes the solution of Ax=B, where A is an arbitrary real nxn matrix.
 * The auxiliary vector piv returns the pivoting row exchange data.
 * */
int LP_gen_solve(double *A, int n, double *B, int nrhs, int *piv)
{
	int info;
	dgesv_(&n, &nrhs, A, &n, piv, B, &n, &info);
	return info;
}




/*********************************************************************************************
logDetFromCholesky - returns the log determinant of an nxn matrix given its cholesky 
                     decomposition m 
*********************************************************************************************/
double logDetFromCholesky(const double *m, int n) {
	double log_det = 0.0;
	int i = 0;
	while (n > 0) {
		log_det += log(m[i]);
		i+= n;
		n--;
	}
	log_det *= 2;
	return (log_det);
}




/*********************************************************************************************
matrix_multiply - stores A%*%B in the matrix ans
*********************************************************************************************/
void matrix_multiply(const double *A, int nrowA, int ncolA, const double *B, int ncolB, double *ans) {
	int i, j, p;
	double total;
	for (i = 0; i < nrowA; i++) {
		for (j = 0; j < ncolB; j++) {
			total = 0.0000;
			for (p = 0; p < ncolA; p++) {
				total += A[i*ncolA+p] * B[p*ncolB + j];
			}
			ans[i*ncolB + j] = total;
		}
	}
}

// multiplies all elements of packed nxn matrix m by c
void packed_matrix_scale_const(double *m, const double c, const int n) {
	int i;
	for (i = 0; i < n*(n+1)/2; i++) {
		m[i] *= c;
	}
}

/***************************************************
static void scale_packed_corMatrix(double *m, const double c, const double *reps, const int n)
{
	int i, j;
	for (i=0; i<n; i++) {
		for (j=0; j<=i; j++) {
		//	printout("m = %f, c = %f, reps[i] = %f,  reps[j] = %f\n", m[i+(j*(2*n-j-1))/2], c, reps[i], reps[j]); 
			m[i+(j*(2*n-j-1))/2] *= (c/ sqrt(reps[i]*reps[j])   );
		}
	}
}
*******************************************/
/***********************************************************************
cbind - implementation of cbind from R; stores the matrix {X1 X2} in A
***********************************************************************/
void cbind(const double *X1, const double *X2, double *A, int ncolX1, int ncolX2, int nrowA) {
	int i;
	int j;
	int ncolA = ncolX1 + ncolX2;

	for (i = 0; i < nrowA; i++) {
		for (j = 0; j < ncolX1; j++) {
			MATRIX_SET(A,i,j,ncolA,  MATRIX_GET(X1,i,j,ncolX1));

		}
		for (j = 0; j< ncolX2; j++) {
			MATRIX_SET(A,i,j+ncolX1,ncolA, MATRIX_GET(X2,i,j,ncolX2));
		}
	}

}



/********************************************************
transpose - stores transpose of matrix A in Aprime 
********************************************************/
void transpose(const double *A, double *Aprime, int nrowA, int ncolA) {
	int i, j;
	for (i = 0; i < nrowA; i++) {
		for (j = 0; j < ncolA; j++) {
			MATRIX_SET(Aprime, j,i,nrowA, MATRIX_GET(A,i,j,ncolA));
		}
	}
}

// calculates v1 - v2 and stores results in v1
void vectorSubtract(double *v1, const double *v2, int n) {
	int i;
	for (i = 0; i < n; i++) {
		v1[i] = v1[i] - v2[i];
	}
}
// copies vector v1 to vector v2
void vectorCopy(const double *v1, double *v2, const int n) {
	int i;
	for (i = 0; i < n; i++) {
		v2[i]=v1[i];
	}
}
// calculates the sample variance of the first n elements in v
double vectorVariance(const double *v, const int n) {
	int i;
	double mean = 0.0;
	for (i = 0; i < n; i++) {
		mean += v[i];
	}
	mean = mean / n;
	
	double sse = 0.0;
	double dev;
	for (i = 0; i < n; i++) {
		dev = v[i] - mean;
		sse += (dev*dev);
	}
	return sse / (double) (n-1);
}

/*****************************************************************
printMatrix
*****************************************************************/
static void printMatrix(const char *fmt, const double *a, int nrow, int ncol, char *header)
{
	if (header != NULL) printout("%s\n", header);
	int i, j;
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncol; j++)
		//	printout(fmt, a[i+n*j]);
			printout(fmt, a[i*ncol+j]);
		//putchar('\n');
	}
}
/*****************************************************************************************
print_packed_matrix - for symmetric matrices where only lower (upper) triangle is stored
******************************************************************************************/
static void printPackedMatrix(const char *fmt, const double *a, int n)
{
	int i, j;
	for (i=0; i<n; i++) {
		for (j=0; j<=i; j++)
			printout(fmt, a[i+(j*(2*n-j-1))/2]);
		//putchar('\n');
	}
}


/*************************************************************************
printVector - prints the first m elements of vector v to 1 row of output
*************************************************************************/
static void printVector(const char *fmt, const double *b, int m)
{
	int i;
	for (i=0; i<m; i++)
		printout(fmt, b[i]);
	//putchar('\n');
}




/***********************************************************************
xprimeA - stores the matrix x'A in the vector ans
	  double *x - pointer to vector of length n (1xn matrix)
	  double *A - pointer to an nxm matrix 
	  double *ans - pointer to vector of length m which stores the answer
	  int nrowA - number of rows of A
	  int ncolsA - number of columns of A
***********************************************************************/
void xprimeA(const double *x, const double *A, double *ans, int nrowA, int ncolA) {
	int i,j;
	double total = 0.0;
	for (j = 0; j < ncolA; j++) {
		total = 0.0;
		for (i = 0; i < nrowA; i++) {
			 total += MATRIX_GET(A,i,j,ncolA) * x[i];
		}
		ans[j] = total;
	}
}


/**************************************************************************
dotprod - calculates and returns the dot product of vectors v1 and v2 
**************************************************************************/
double dotprod(const double *v1, const double *v2, int n) {
	int i;
	double total = 0.0;
	for (i = 0; i < n; i++) {
		total += v1[i]*v2[i];
	}
	return total;
}

#endif



