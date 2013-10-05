#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "i_malloc.h"

void init_allocator(void* malloc_fn, void* calloc_fn, void* realloc_fn, void* free_fn) {
  i_malloc = malloc_fn;
  i_calloc = calloc_fn;
  i_realloc = realloc_fn;
  i_free = free_fn; 
}

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


/* allocation and de-allocation of aligned data for use as matrix storage */
double* alloc_doubles(unsigned long n) {
  return (double*)MKL_malloc(n*sizeof(double), 64);
}

void free_doubles(double* ptr) {
  MKL_free(ptr);
}


/* matrix copy operations */

void copy_matrix(const double* from_data, unsigned long num_rows, unsigned long num_cols, double* to_data) {
  char uplo = 'A'; /* copy whole matrix */
  MKL_INT rows = (MKL_INT)num_rows;
  MKL_INT cols = (MKL_INT)num_cols;
  /*void dlacpy( const char* uplo, const MKL_INT* m, const MKL_INT* n,
                 const double* a, const MKL_INT* lda, double* b, const MKL_INT* ldb );*/
  dlacpy(&uplo, &rows, &cols,
         from_data, &rows, to_data, &rows);
}

void copy_upper_triangle(const double* from_data, unsigned long num_rows, unsigned long num_cols, double* to_data) {
  char uplo = 'U'; /* copy upper triangle only */
  MKL_INT rows = (MKL_INT)num_rows;
  MKL_INT cols = (MKL_INT)num_cols;
  /*void dlacpy( const char* uplo, const MKL_INT* m, const MKL_INT* n,
                 const double* a, const MKL_INT* lda, double* b, const MKL_INT* ldb );*/
  dlacpy(&uplo, &rows, &cols,
         from_data, &rows, to_data, &rows);
}


