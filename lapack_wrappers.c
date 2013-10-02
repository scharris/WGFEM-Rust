#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

void printm(double* a, lapack_int m, lapack_int n, lapack_int lda) {
	lapack_int i, j;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}

void printiv(lapack_int* a, lapack_int n) {
	lapack_int j;
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}


lapack_int solve_symmetric_as_col_maj_with_ut_sys(double* a,
                                                  lapack_int n,
                                                  double* b,
                                                  lapack_int nrhs,
                                                  lapack_int* ipiv) {
  lapack_int info;
  
  /* printm(a, n, n, n); printm(b, n, nrhs, n); */

  info = LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', n, nrhs, a, n, ipiv, b, n);

  /* printm(b, n, nrhs, n); printm(a, n, n, n); printiv(ipiv, n); */
} 

