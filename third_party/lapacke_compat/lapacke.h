// Minimal LAPACKE declarations needed by vendored Algoim.
//
// Some conda-forge LAPACK builds export LAPACKE symbols from liblapack but do
// not ship lapacke.h. Algoim only calls the two routines declared here.

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifndef LAPACK_ROW_MAJOR
#define LAPACK_ROW_MAJOR 101
#endif

typedef int lapack_int;

lapack_int LAPACKE_dgesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, double* a,
                          lapack_int lda, double* s, double* u,
                          lapack_int ldu, double* vt, lapack_int ldvt,
                          double* superb);

lapack_int LAPACKE_dggevx(int matrix_layout, char balanc, char jobvl,
                          char jobvr, char sense, lapack_int n, double* a,
                          lapack_int lda, double* b, lapack_int ldb,
                          double* alphar, double* alphai, double* beta,
                          double* vl, lapack_int ldvl, double* vr,
                          lapack_int ldvr, lapack_int* ilo,
                          lapack_int* ihi, double* lscale, double* rscale,
                          double* abnrm, double* bbnrm, double* rconde,
                          double* rcondv);

#ifdef __cplusplus
}
#endif
