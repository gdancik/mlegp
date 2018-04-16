#ifndef H_LAPACK_HEADERS_H
#define H_LAPACK_HEADERS_H

extern void dgesv_(int *n, int *nrhs, double *A, int *lda, int *piv,
		double *B, int *ldb, int *info);

extern void dppsv_(char *uplo, int *n, int *nrhs, double *A, 
		double *B, int *ldb, int *info);

extern void dsyevr_(char *jobz, char *range, char *uplo,
		int *n, double *a, int *lda,
		double *vl, double *vu, int *il, int *iu,
		double *abstol, int *m, double *w,
		double *z, int *ldz, int *isuppz,
		double *work, int *lwork, int *iwork, int *liwork,
		int *info);

#endif /* H_LAPACK_HEADERS_H */
