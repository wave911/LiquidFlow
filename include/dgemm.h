#ifdef __cplusplus
extern "C" {
#endif

#include "f2c.h"

int dgemm_(char *TransA, char *TransB, integer *m, integer *n, integer *k, doublereal *alpha, doublereal *A, integer *lda, doublereal *B, integer *ldb, doublereal *beta, doublereal *C, integer *ldc);

#ifdef __cplusplus
}
#endif
