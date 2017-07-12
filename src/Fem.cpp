#include "Fem.h"
#include <iostream>
#ifdef LIBBLASLAPACK
#include "f2c.h"
#include "dgemv.h"
#include "dgesv.h"
#endif
#ifdef APPLE
#include <Accelerate/Accelerate.h>
#endif

using namespace std;

void CFem::dgemv(char *trans, int m, int n, real_t alpha, real_t *a, 
				   int lda, real_t *x, int incx, real_t beta, real_t *y, int incy) {
#ifdef LIBBLASLAPACK
	dgemv_(trans, (integer*)&m, (integer*)&n, (doublereal*)&alpha,
				(doublereal*)a, (integer*)&lda, (doublereal*)x,
				(integer*)&incx, (doublereal*)&beta, (doublereal*)y,
				(integer*)&incy);
#endif
#ifdef APPLE
	cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
			    m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif
}

void CFem::dgesv(int n, real_t *M, real_t *B ) {
	int status = 0;
	int NRHS = 1;
	int* IPIV = new int[n];
	int N = n;
    int lda = n;
    int ldb = n;
#ifdef LIBBLASLAPACK
	dgesv_((integer*)&N, (integer*)&NRHS, M, (integer*)&lda, (integer*)IPIV, B, (integer*)&ldb, (integer*)&status);
#endif
#ifdef APPLE
	dgesv_((__CLPK_integer*)&N, (__CLPK_integer*)&NRHS, M, (__CLPK_integer*)&lda, (__CLPK_integer*)IPIV, B, (__CLPK_integer*)&ldb, (__CLPK_integer*)&status);
#endif
	delete [] IPIV;
}




