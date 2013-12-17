#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
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


/* allocation and de-allocation of aligned data for use as matrix storage */
double* alloc_doubles(unsigned long n) {
  return (double*)MKL_malloc(n*sizeof(double), 64);
}

void free_doubles(double* ptr) {
  MKL_free(ptr);
}

/* allocation and de-allocation of aligned data for use as matrix storage */
lapack_int* alloc_ints(unsigned long n) {
  return (lapack_int*)MKL_malloc(n*sizeof(lapack_int), 64);
}

void free_ints(lapack_int* ptr) {
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

/* Dense matrix system solver. */
lapack_int solve_symmetric_as_col_maj_with_ut_sys
           (double* a, lapack_int n, double* b, lapack_int nrhs, lapack_int* ipiv) {
  return LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', n, nrhs, a, n, ipiv, b, n);
} 


/* Sparse symmetric matrix system solver. */
MKL_INT solve_sparse_symmetric_as_ut_csr3(MKL_INT n, const MKL_INT* ia, const MKL_INT* ja, const double* a,
                                          const double* b, MKL_INT nrhs,
                                          double* x,
                                          unsigned num_cpu_cores) {

  MKL_INT mtype = -2; /* symmetric indefinite */
  void *pt[64];
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;

  MKL_INT i, i_un;
  double d_un; /* "*_un" for unused params */

  for (i = 0; i<64; i++) { iparm[i] = 0; }
  iparm[0] = 1;  /* Not all defaults */
  iparm[1] = 2;  /* Fill-in reordering from METIS */
  iparm[7] = 15; /* Max numbers of iterative refinement steps. 0 also means 2 iterations but does not seem to allow for early stopping.  */
  iparm[9] = 8;  /* Perturb the pivot elements with 1E-8 which is the default for symmetric matrices. */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[12] = 1; /* Maximum weighted matching algorithm (default off for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
  iparm[20] = 1; /* pivoting method, Bunch-Kaufman is recommended for symmetric indefinite matrices */
  iparm[23] = num_cpu_cores > 8 ? 1 : 0; /* Use two level parallel factorization algorithm. */
  iparm[26] = 1; /* Check matrix. TODO: Unset after testing. */
  iparm[34] = 1; /* Use 0-based row and column numbers within ia and ja arrays. */
 
  maxfct = 1;    /* Leave this at 1. Number of numerical factorizations to keep in memory */
  mnum = 1;      /* Leave this at 1. Which factorization of the above to use in the solving step. */
  msglvl = 0;    /* No statistical information output. */
  error = 0;
  
  /* Required initialization for internal data pointer. */
  for (i = 0; i<64; i++) { pt[i] = 0; }

  /* Reordering and Symbolic Factorization. */
  phase = 11;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &i_un, &nrhs, iparm, &msglvl, &d_un, &d_un, &error);
  
  if (error != 0) { fprintf (stderr, "\nERROR during symbolic factorization: %d", error); return error; }
  
  /* Numerical factorization. */
  phase = 22;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &i_un, &nrhs, iparm, &msglvl, &d_un, &d_un, &error);
  
  if (error != 0) { fprintf (stderr, "\nERROR during numerical factorization: %d", error); return error; }

  /* Back substitution and iterative refinement. */
  phase = 33;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &i_un, &nrhs, iparm, &msglvl, b, x, &error);

  if (error != 0) { fprintf (stderr, "\nERROR during solution: %d", error); return error; }

  /* Release resources. */
  phase = -1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &d_un, ia, ja, &i_un, &nrhs, iparm, &msglvl, &d_un, &d_un, &error);
  return 0;
}

/* Sparse structurally symmetric matrix system solver. */
MKL_INT solve_sparse_structurally_symmetric_csr3(MKL_INT n, const MKL_INT* ia, const MKL_INT* ja, const double* a,
                                                 const double* b, MKL_INT nrhs,
                                                 double* x,
                                                 unsigned num_cpu_cores) {

  MKL_INT mtype = 1; /*  1 == Real structurally symmetric matrix. */
  void *pt[64];
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  
  MKL_INT i, i_un;
  double d_un;  /* "*_un" for unused params */

  for (i=0; i<64; i++) { iparm[i] = 0; }
  iparm[0] = 1;  /* Not all defaults */
  iparm[1] = 2;  /* Fill-in reordering from METIS */
  iparm[7] = 20;  /* Max numbers of iterative refinement steps. 0 also means 2 iterations but does not seem to allow for early stopping.  */
  iparm[9] = 13; /* Perturb the pivot elements with 1E-13 which is the default for nonsymmetric matrices. */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS. */
  iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for nonsymmetric matrices). */
  iparm[23] = num_cpu_cores > 8 ? 1 : 0; /* Use two level parallel factorization algorithm. */
  iparm[26] = 1; /* Check matrix. TODO: Unset after testing. */
  iparm[34] = 1; /* Use 0-based row and column numbers within ia and ja arrays. */

  maxfct = 1;    /* Leave this at 1. Number of numerical factorizations to keep in memory. */
  mnum = 1;      /* Leave this at 1. Which factorization of the above to use in the solving step. */
  msglvl = 1;    /* No statistical information output. */
  error = 0;
  
  /* Required initialization for MKL internal use data. */
  for (i=0; i<64; i++) { pt[i] = 0; }

  /*  Reordering and Symbolic Factorization. */
  phase = 11;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &i_un, &nrhs, iparm, &msglvl, &d_un, &d_un, &error);
  
  if (error != 0) { fprintf (stderr, "\nERROR during symbolic factorization: %d", error); return error; }

  /* Numerical factorization. */
  phase = 22;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &i_un, &nrhs, iparm, &msglvl, &d_un, &d_un, &error);
 
  if (error != 0) { fprintf ("\nERROR during numerical factorization: %d", error); return error; }

  /* Back substitution and iterative refinement. */
  phase = 33;

  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &i_un, &nrhs, iparm, &msglvl, b, x, &error);
 
  if (error != 0) { fprintf ("\nERROR during solution: %d", error); return error; }

  /* Release resources. */
  phase = -1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &d_un, ia, ja, &i_un, &nrhs, iparm, &msglvl, &d_un, &d_un, &error);
  return 0;
}

