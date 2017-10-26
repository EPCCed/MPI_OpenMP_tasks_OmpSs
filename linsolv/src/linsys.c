/*
 * linsys.c
 */

#include "linsys.h"

#include <math.h>
#include "matrix.h"
#include "util.h"

static char *linsolv_names[NBSLS_SOLVERS] =
{
  "Block point implicit",
  "Block Jacobi",
  "Block Gauss-Seidel",
  "Block symmetric Gauss-Seidel"
};

/*******************************************************************************
* Swaps two rows of a matrix.
*
* @param[in,out] matrix The matrix.
* @param[in] rowOne Index of first row which is swapped with second row.
* @param[in] rowTwo Index of second row which is swapped with first row.
* @warning rowOne and rowTwo must be valid row indices of matrix.
*******************************************************************************/
static void swap_rows(Matrix matrix, int rowOne, int rowTwo);

/*******************************************************************************
* Solves linear system Lx=b where L is a lower triangular quadratic matrix
* and has entry 1 on the diagonal, i.e. L[i][i] = 1.0, i = 0,...,n-1.
*
* @param[in] L     Lower triangular quadratic matrix L, L must have entry 1 on
*                  the diagonal.
* @param[in,out] b Right hand side vector b, solution of Lx=b.
* @return 0
*******************************************************************************/
static int substitute_forward(const Matrix L, double *b);

/*******************************************************************************
* Solves linear system Ux=b where U is upper triangular quadratic matrix.
*
* @param[in] U     Upper triangular quadratic matrix U, U must be invertible.
* @param[in,out] b Right hand side vector b, solution of Ux=b.
* @return 0
*******************************************************************************/
static int substitute_backward(const Matrix U, double *b);

/*******************************************************************************
*
*******************************************************************************/
void delete_bsm_linear_system(BSMLinSys bsmls)
{
  int row;

  for(row = 0; row < bsmls.A.num_rows; row++)
    check_free(bsmls.swap[row]);
  check_free(bsmls.swap);

  deleteMatrix(bsmls.rhs);

  delete_block_sparse_matrix(bsmls.A);
} /* delete_bsm_linear_system() */

/*******************************************************************************
*
*******************************************************************************/
char *get_linsolv_name(BSLinSysMethod solverID)
{
  return linsolv_names[solverID];
} /* get_linsolv_name() */

/*******************************************************************************
*
*******************************************************************************/
void lu_decomp(Matrix A, int *swap)
{
  const int rows = A.rows;
  int maxPos, i, j, k;

  for(k = 0; k < rows - 1; k++)
  {
    maxPos = k;

    for(i = k+1; i < A.rows; i++) /* find max in column k */
      if(fabs(A.m[i][k]) > fabs(A.m[maxPos][k]))
        maxPos = i;

    swap[k] = maxPos;
    swap_rows(A, maxPos, k);

    for(i = k+1; i < rows; i++)
    {
      A.m[i][k] = A.m[i][k] / A.m[k][k];

      for(j = k+1; j < rows; j++)
        A.m[i][j] = A.m[i][j] - A.m[i][k] * A.m[k][j];
    }
  }
} /* lu_decomp() */

/*******************************************************************************
*
*******************************************************************************/
int lu_solve(const Matrix A, double *b, const int *swap)
{
  int i;
  double help;

  for(i = 0; i < A.rows - 1; i++)
  {
    help = b[i];
    b[i] = b[swap[i]];
    b[swap[i]] = help;
  }

  if(substitute_forward(A, b) == 0)
    if(substitute_backward(A, b) == 0)
      return 0;
  return 1;
} /* lu_solve() */

/*******************************************************************************
*
*******************************************************************************/
static void swap_rows(Matrix matrix, int rowOne, int rowTwo)
{
  double *vector = matrix.m[rowOne];

  matrix.m[rowOne] = matrix.m[rowTwo];
  matrix.m[rowTwo] = vector;
} /* swap_rows() */

/*******************************************************************************
*
*******************************************************************************/
static int substitute_forward(const Matrix L, double *b)
{
  int i, j;
  double help;

  for(i = 0; i < L.rows; i++)
  {
    help = 0.0;
    for(j = 0; j < i; j++)
      help += L.m[i][j] * b[j];
    b[i] = b[i] - help;
  }
  return 0;
} /* substitute_forward() */

/*******************************************************************************
*
*******************************************************************************/
static int substitute_backward(const Matrix U, double *b)
{
  int i, j;
  double help;

  for(i = U.rows - 1; i >= 0; i--)
  {
    help = 0.0;
    for(j = U.rows - 1; j > i; j--)
      help += U.m[i][j] * b[j];
    b[i] = (b[i] - help) / U.m[i][i];
  }
  return 0;
} /* substitute_backward() */
