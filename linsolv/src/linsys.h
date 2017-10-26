/*
 * linsys.h
 */

#ifndef SRC_LINSYS_H_
#define SRC_LINSYS_H_

#include "matrix.h"

/*******************************************************************************
* Constants of possible decomposition states of matrices within block diagonal
* linear systems, block tridiagonal linear systems and linear systems having a
* block sparse matrix.
*******************************************************************************/
typedef enum
{
  /* Matrix is not decomposed */
  NOT_DECOMPOSED = 0,
  /* The diagonal block of a block sparse matrix are decomposed by
     LU with pivoting */
  BSM_DIAG_LU_PIVOT,
  /* In case line information is available the corresponding tridiagonal
     blocks of a block sparse matrix are decomposed by block LU and diagonal
     elements not corresponding to a line are decomposed by LU with pivoting */
  BSM_TRI_DIAG_LU_PIVOT
} MatDecompStat;

/*******************************************************************************
* Constants of possible iterative solution methods for solving block sparse
* linear systems.
*******************************************************************************/
typedef enum
{
  BSLS_POINT_IMPLICIT = 0,
  BSLS_JACOBI,
  BSLS_GAUSS_SEIDEL,
  BSLS_SYMM_GAUSS_SEIDEL,
  NBSLS_SOLVERS
} BSLinSysMethod;

/*******************************************************************************
* Structure defining a block sparse matrix linear system. The block sparse
* matrix is given in Block Compressed Sparse Row (BCSR) format, the right
* hand side is given by the Matrix rhs.
*******************************************************************************/
typedef struct
{
  /* Block sparse matrix of linear system */
  BlockSparseMatrix A;

  /* Right hand side of linear system */
  Matrix rhs;

  /* Memory for LU-Pivot decomposition of diagonal elements */
  int **swap;

  /* Status of diagonal blocks of BlockSparseMatrix A */
  MatDecompStat decomposition_state;

  /* Diagonal entry is stored first for each column (1) or not (0) */
  int diag_first;
} BSMLinSys;

/*******************************************************************************
* Frees memory which was allocated for a block sparse matrix linear system in
* Block Compressed Sparse Row (BCSR) format.
*
* @param[in] bsmls Block sparse matrix linear system for which memory is freed.
*******************************************************************************/
void delete_bsm_linear_system(BSMLinSys bsmls);

/*******************************************************************************
* Get (printable) name for the given linear solution method.
*******************************************************************************/
char *get_linsolv_name(BSLinSysMethod solverID);

/*******************************************************************************
* Compute by Gauss-Algorithm with column-pivoting an LU-decomposition of
* quadratic matrix A with a lower triangular matrix L with 1's on the
* diagonal and an upper triangular matrix U such that PA = LU where P is
* the permutation matrix arising from pivoting.
*
* @param[in,out] A    Given quadratic matrix A,  Matrix A is overwritten
*                     in the lower triangular part by the matrix L and in the
*                     upper triangular part by the matrix U.
*                     To get L one has to add the 1's on the diagonal
* @param[in,out] swap Array representing the permutation matrix P. Array swap
*                     has dimension A.rows, i.e. swap[A.rows],
*                     out: swap represents matrix P, that is swap[0] is the row
*                     which needs to be swapped with row 0 of A, swap[1] is the
*                     row which ...
* @warning Gauss-Algorithm with column pivoting must be applicable to
*          matrix A.
*******************************************************************************/
void lu_decomp(Matrix A, int *swap);

/*******************************************************************************
* Solves a given linear system Ax = b with the help of a LU-decomposition with
* column pivoting and forward and backward substitution. It is assumed that
* the matrix A is already decomposed, e.g. by the function lu_decomp().
*
* @param[in] A     Given quadratic matrix A which is pivot-LU decomposed.
* @param[in,out] b Given right hand side b of linear system, out: b has
*                  solution of linear system Ax=b.
* @param[in] swap  Array describing the permutation due to pivoting.
* @return 0 Linear system was successfully solved.
* @return 1 Linear system could not be solved.
* @warning Gauss-Algorithm with column pivoting must be applicable to matrix.
* @warning The array swap is of dimension A.rows and must represent the
*          permutation matrix, e.g. computed by lu_decomp().
*******************************************************************************/
int lu_solve(const Matrix A, double *b, const int *swap);

#endif /* SRC_LINSYS_H_ */
