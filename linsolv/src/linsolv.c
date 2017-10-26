/*
 * linsolv.c
 */

#include "linsolv.h"

#include <omp.h>
#include "setup_comm.h"
#include "exchange_matrix.h"
#include "linsys.h"
#include "matrix.h"
#include "util.h"

#if USE_TASKLOOP == 0
static const int chunksize = SET_CHUNKSIZE;
#endif

/*******************************************************************************
* Solves the block sparse linear system A*x0 = rhs only and directly for the
* decoupled small scale linear systems corresponding to the diagonal
* entries of A.
* Here A denotes the given block sparse matrix, rhs is the right hand
* side of the linear system and x0 is overwritten by the solution of the
* corresponding decoupled linear systems.
*
* @param[in] linsys Given linear system with LU decomposition of the diagonal
*                   entries of the block sparse matrix.
* @param[in,out] x0 On input, initial guess required for iteration.
*                   On output, solution of given linear system corresponding to
*                   the diagonal entries.
*******************************************************************************/
static void point_implicit(BSMLinSys *linsys, Matrix x0);

/*******************************************************************************
* Solves approximately the block sparse linear system A*x0 = rhs by a Jacobi
* method. Here A denotes the given block sparse matrix, rhs is the right hand
* side of the linear system and x0 is the initial guess. The maximum number
* of iterations may be defined and needs to be larger than zero.
*
* @param[in] linsys   Given linear system with LU decomposition of the diagonal
*                     entries of the block sparse matrix.
* @param[in] parallel_data Data structure required for parallelization.
* @param[in,out] x0   On input, initial guess required for iteration.
*                     On output, approximate solution of given linear system.
* @param[in] max_iter Maximum number of iterations.
*******************************************************************************/
static void jacobi(BSMLinSys *linsys, CommMap *parallel_data, Matrix x0,
                   const int max_iter);

/*******************************************************************************
* Solves approximately the block sparse linear system A*x0 = rhs by a
* Gauss-Seidel method. Here A denotes the given block sparse matrix, rhs is
* the right hand side of the linear system and x0 is the initial guess. The
* maximum number of iterations may be defined and needs to be larger than zero.
*
* @param[in] linsys   Given linear system with LU decomposition of the diagonal
*                     entries of the block sparse matrix.
* @param[in] parallel_data data structure required for parallelization
* @param[in,out] x0   On input, initial guess required for iteration.
*                     On output, approximate solution of given linear system.
* @param[in] max_iter Maximum number of iterations.
*******************************************************************************/
static void gauss_seidel(BSMLinSys *linsys, CommMap *parallel_data, Matrix x0,
                         int max_iter);

/*******************************************************************************
* Solves approximately the block sparse linear system A*x0 = rhs by a symmetric
* Gauss-Seidel method. Here A denotes the given block sparse matrix, rhs is
* the right hand side of the linear system and x0 is the initial guess. The
* maximum number of iterations may be defined and needs to be larger than zero.
*
* @param[in] linsys   Given linear system with LU decomposition of the diagonal
*                     entries of the block sparse matrix.
* @param[in] parallel_data Data structure required for parallelization.
* @param[in,out] x0   On input, initial guess required for iteration.
*                     On output, approximate solution of given linear system.
* @param[in] max_iter Maximum number of iterations.
*******************************************************************************/
static void symm_gauss_seidel(BSMLinSys *linsys, CommMap *parallel_data,
                              Matrix x0, int max_iter);

#if USE_TASKLOOP == 0
/*******************************************************************************
*
*******************************************************************************/
void linsys_lu_decomp_diag(BSMLinSys *bsmls)
{
  BlockSparseMatrix A = bsmls->A;
  int **swap = bsmls->swap;
  int *row_ptr = A.row_ptr;
  Matrix *entries = A.entries;
  const int num_rows = A.num_rows;

  CHECK(bsmls->decomposition_state == NOT_DECOMPOSED);
  CHECK(bsmls->diag_first);

# pragma omp parallel default(none) shared(entries, row_ptr, swap)
  {
#   pragma omp single
    {
      for(int i = 0; i < num_rows; i += chunksize)
      {
#       pragma omp task
        {
          for(int row = i; row < MIN(i + chunksize, num_rows); row++)
            lu_decomp(entries[row_ptr[row]], swap[row]);
        }
      }
    } /* end single */
  } /* end of parallel region */

  bsmls->decomposition_state = BSM_DIAG_LU_PIVOT;
} /* linsys_lu_decomp_diag() */
#endif /* USE_TASKLOOP == 0 */

#if USE_TASKLOOP == 1
/*******************************************************************************
*
*******************************************************************************/
void linsys_lu_decomp_diag(BSMLinSys *bsmls)
{
  BlockSparseMatrix A = bsmls->A;
  int **swap = bsmls->swap;
  int *row_ptr = A.row_ptr;
  Matrix *entries = A.entries;
  const int num_rows = A.num_rows;

  CHECK(bsmls->decomposition_state == NOT_DECOMPOSED);
  CHECK(bsmls->diag_first);

# pragma omp parallel default(none) shared(entries, row_ptr, swap)
  {
#   pragma omp single
    {
#     pragma omp taskloop GRAINSIZECLAUSE
      for(int row = 0; row < num_rows; row++)
        lu_decomp(entries[row_ptr[row]], swap[row]);
    } /* end single */
  } /* end of parallel region */

  bsmls->decomposition_state = BSM_DIAG_LU_PIVOT;
} /* linsys_lu_decomp_diag() */
#endif /* USE_TASKLOOP == 1 */

/*******************************************************************************
*
*******************************************************************************/
void linsolv(const BSLinSysMethod linear_solver, const int num_sweeps,
             BSMLinSys *linsys, Matrix approximate_solution, CommMap *comap)
{
  /*----------------------------------------------------------------------------
  | Depending on the choice of the iterative linear solution methodology, one
  | of the following list is chosen. All these are some kind of
  | block Jacobi or block Gauss-Seidel method.
  ----------------------------------------------------------------------------*/
  switch(linear_solver)
  {
    case BSLS_POINT_IMPLICIT:
    {
      point_implicit(linsys, approximate_solution);
      break;
    }
    case BSLS_JACOBI:
    {
      jacobi(linsys, comap, approximate_solution, num_sweeps);
      break;
    }
    case BSLS_GAUSS_SEIDEL:
    {
      gauss_seidel(linsys, comap, approximate_solution, num_sweeps);
      break;
    }
    case BSLS_SYMM_GAUSS_SEIDEL:
    {
      symm_gauss_seidel(linsys, comap, approximate_solution, num_sweeps);
      break;
    }
    default:
    {
      break;
    }
  } /* switch(linear_solver) */
} /* linsolv() */

#if USE_TASKLOOP == 0
/*******************************************************************************
*
*******************************************************************************/
static void point_implicit(BSMLinSys *linsys, Matrix x0)
{
  BlockSparseMatrix A = linsys->A;
  int *row_ptr = A.row_ptr;
  Matrix *entries = A.entries;
  const int num_rows = A.num_rows;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

  /*----------------------------------------------------------------------------
  | only loop over linear systems corresponding to the diagonal entries
  ----------------------------------------------------------------------------*/
# pragma omp parallel default(none) shared(entries, row_ptr, rhs, swap, x0)
  {
#   pragma omp single
    {
      for(int i = 0; i < num_rows; i += chunksize)
      {
#       pragma omp task
        {
          for(int row = i; row < MIN(i + chunksize, num_rows); row++)
          {
            for(int j = 0; j < rhs.cols; j++)
              x0.m[row][j] = rhs.m[row][j];

            lu_solve(entries[row_ptr[row]], x0.m[row], swap[row]);
          }
        }
      }
    } /* end single */
  } /* end of parallel region */
} /* point_implicit() */
#endif /* USE_TASKLOOP == 0 */

#if USE_TASKLOOP == 1
/*******************************************************************************
*
*******************************************************************************/
static void point_implicit(BSMLinSys *linsys, Matrix x0)
{
  BlockSparseMatrix A = linsys->A;
  int *row_ptr = A.row_ptr;
  Matrix *entries = A.entries;
  const int num_rows = A.num_rows;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

  /*----------------------------------------------------------------------------
  | only loop over linear systems corresponding to the diagonal entries
  ----------------------------------------------------------------------------*/
# pragma omp parallel default(none) shared(entries, row_ptr, rhs, swap, x0)
  {
#   pragma omp single
    {
#     pragma omp taskloop GRAINSIZECLAUSE
      for(int row = 0; row < num_rows; row++)
      {
        for(int j = 0; j < rhs.cols; j++)
          x0.m[row][j] = rhs.m[row][j];

        lu_solve(entries[row_ptr[row]], x0.m[row], swap[row]);
      }
    } /* end single */
  } /* end of parallel region */
} /* point_implicit() */
#endif /* USE_TASKLOOP == 1 */

#if USE_TASKLOOP == 0
/*******************************************************************************
*
*******************************************************************************/
static void jacobi(BSMLinSys *linsys, CommMap *parallel_data, Matrix x0,
                   const int max_iter)
{
  BlockSparseMatrix A = linsys->A;
  int *row_ptr = A.row_ptr;
  Matrix *entries = A.entries;
  int *col_idx = A.col_index;
  const int num_rows = A.num_rows;
  const int brows = A.block_size_row;
  const int bcols = A.block_size_col;
  Matrix rhs = linsys->rhs;
  const int rhs_rows = rhs.rows;
  int **swap = linsys->swap;

  Matrix help  = generateMatrix(rhs.rows, rhs.cols);

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);
# pragma omp parallel default(none) shared(row_ptr, col_idx, entries, swap, \
                                           rhs, x0, help, parallel_data)
  {
#   pragma omp single
    {
      int row, idx, col, i, j, k;
      int iter = 0;

      /*------------------------------------------------------------------------
      | perform a Jacobi iteration
      ------------------------------------------------------------------------*/
      do
      {
        iter++;

        /* copy old solution */
#       pragma omp taskgroup
        {
          for(k = 0; k < rhs_rows; k += chunksize)
          {
#           pragma omp task
            {
              for(i = k; i < MIN(k + chunksize, rhs_rows); i++)
              {
                for(j = 0; j < rhs.cols; j++)
                  help.m[i][j] = x0.m[i][j];
              }
            }
          }
        } /* end taskgroup */

        /* calculate new solution */
#       pragma omp taskgroup
        {
          for(k = 0; k < num_rows; k += chunksize)
          {
#           pragma omp task
            {
              for(row = k; row < MIN(k + chunksize, num_rows); row++)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  col = col_idx[idx];
                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= entries[idx].m[i][j] * help.m[col][j];
                }

                lu_solve(entries[row_ptr[row]], x0.m[row], swap[row]);
              }
            } /* end task */
          } /* end for */
        } /* end taskgroup */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

      } while(iter < max_iter);
    } /* end single */
  } /* end of parallel region */

  deleteMatrix(help);
} /* jacobi() */
#endif /* USE_TASKLOOP == 0 */

#if USE_TASKLOOP == 1
/*******************************************************************************
*
*******************************************************************************/
static void jacobi(BSMLinSys *linsys, CommMap *parallel_data, Matrix x0,
                   const int max_iter)
{
  BlockSparseMatrix A = linsys->A;
  int *row_ptr = A.row_ptr;
  Matrix *entries = A.entries;
  int *col_idx = A.col_index;
  const int num_rows = A.num_rows;
  const int brows = A.block_size_row;
  const int bcols = A.block_size_col;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  Matrix help  = generateMatrix(rhs.rows, rhs.cols);

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

# pragma omp parallel default(none) shared(row_ptr, col_idx, entries, swap, \
                                           rhs, x0, help, parallel_data)
  {
#   pragma omp single
    {
      int row, idx, col, i, j;
      int iter = 0;

      /*--------------------------------------------------------------------------
      | perform a Jacobi iteration
      --------------------------------------------------------------------------*/
      do
      {
        iter++;

        /* copy old solution */
#       pragma omp taskloop GRAINSIZECLAUSE
        for(i = 0; i < rhs.rows; i++)
        {
          for(j = 0; j < rhs.cols; j++)
            help.m[i][j] = x0.m[i][j];
        }

        /* calculate new solution */
#       pragma omp taskloop GRAINSIZECLAUSE
        for(row = 0; row < num_rows; row++)
        {
          for(i = 0; i < rhs.cols; i++)
            x0.m[row][i] = rhs.m[row][i];

          for(idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
          {
            col = col_idx[idx];
            for(i = 0; i < brows; i++)
              for(j = 0; j < bcols; j++)
                x0.m[row][i] -= entries[idx].m[i][j] * help.m[col][j];
          }

          lu_solve(entries[row_ptr[row]], x0.m[row], swap[row]);
        } /* end taskloop */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

      } while(iter < max_iter);
    } /* end single */
  } /* end of parallel region */

  deleteMatrix(help);
} /* jacobi() */
#endif /* USE_TASKLOOP == 1 */

#if USE_TASKLOOP == 0
/*******************************************************************************
*
*******************************************************************************/
static void gauss_seidel(BSMLinSys *linsys, CommMap *parallel_data, Matrix x0,
                         int max_iter)
{
  BlockSparseMatrix A = linsys->A;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  int *row_ptr = A.row_ptr;
  int *col_idx = A.col_index;
  const int brows = A.block_size_row;
  const int bcols = A.block_size_col;
  const int num_rows = A.num_rows;
  const int rhs_rows = rhs.rows;

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

  Matrix help  = generateMatrix(rhs.rows, rhs.cols);

# pragma omp parallel default(none) shared(A, col_idx, row_ptr, rhs, help, \
                                           parallel_data, swap, max_iter, x0)
  {
    int nthreads = omp_get_num_threads();
    int iter = 0;
    int i, j, k;
    CHECK(nthreads > 0);

    /* for loops with one chunk per thread */
    const int chunkmin = num_rows / nthreads;
    const int rest = num_rows % nthreads;

#   pragma omp single
    {
      /*------------------------------------------------------------------------
      | perform a Gauss-Seidel iteration
      ------------------------------------------------------------------------*/
      do
      {
        iter++;

        /* copy old solution */
#       pragma omp taskgroup
        {
          for(k = 0; k < rhs_rows; k += chunksize)
          {
#           pragma omp task
            {
              for(i = k; i < MIN(k + chunksize, rhs_rows); i++)
              {
                for(j = 0; j < rhs.cols; j++)
                  help.m[i][j] = x0.m[i][j];
              }
            }
          }
        } /* copy */

        /* calculate new solution */
#       pragma omp taskgroup
        {
          for(int tid = 0; tid < nthreads; tid++)
          {
            int mychunksize = chunkmin + (tid < rest);
            int start = tid * chunkmin + MIN(tid, rest);
            int stop = start + mychunksize;

#           pragma omp task
            {
              for(int row = start; row < stop; row++)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(int idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  double *x0_row = NULL;
                  int col = col_idx[idx];

                  if(col >= start && col < stop)
                    x0_row = x0.m[col]; /* use own new solution */
                  else
                    x0_row = help.m[col]; /* use old solution */

                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= A.entries[idx].m[i][j] * x0_row[j];
                } /* end idx loop */

                lu_solve(A.entries[row_ptr[row]], x0.m[row], swap[row]);
              } /* end row loop */
            } /* end task */
          } /* end tid loop */
        } /* calc */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

      } while(iter < max_iter);
    } /* end single */
  } /* end of parallel region */
} /* gauss_seidel() */
#endif /* USE_TASKLOOP == 0 */

#if USE_TASKLOOP == 1
/*******************************************************************************
*
*******************************************************************************/
static void gauss_seidel(BSMLinSys *linsys, CommMap *parallel_data, Matrix x0,
                         int max_iter)
{
  BlockSparseMatrix A = linsys->A;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  int *row_ptr = A.row_ptr;
  int *col_idx = A.col_index;
  const int brows = A.block_size_row;
  const int bcols = A.block_size_col;
  const int num_rows = A.num_rows;
  const int rhs_rows = rhs.rows;

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

  Matrix help  = generateMatrix(rhs.rows, rhs.cols);

# pragma omp parallel default(none) shared(A, col_idx, row_ptr, rhs, help, \
                                           parallel_data, swap, max_iter, x0)
  {
    int nthreads = omp_get_num_threads();
    int iter = 0;
    int i, j;

    CHECK(nthreads > 0);

    /* one chunk per thread */
    const int chunkmin = num_rows / nthreads;
    const int rest = num_rows % nthreads;

#   pragma omp single
    {
      /*------------------------------------------------------------------------
      | perform a Gauss-Seidel iteration
      ------------------------------------------------------------------------*/
      do
      {
        iter++;

        /* copy old solution */
#       pragma omp taskloop GRAINSIZECLAUSE
        for(int row = 0; row < rhs_rows; row++)
        {
          for(j = 0; j < rhs.cols; j++)
            help.m[row][j] = x0.m[row][j];
        }

        /* calculate new solution */
#       pragma omp taskgroup
        {
          for(int tid = 0; tid < nthreads; tid++)
          {
            int mychunksize = chunkmin + (tid < rest);
            int start = tid * chunkmin + MIN(tid, rest);
            int stop = start + mychunksize;

#           pragma omp task
            {
              for(int row = start; row < stop; row++)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(int idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  double *x0_row = NULL;
                  int col = col_idx[idx];

                  if(col >= start && col < stop)
                    x0_row = x0.m[col]; /* use own new solution */
                  else
                    x0_row = help.m[col]; /* use old solution */

                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= A.entries[idx].m[i][j] * x0_row[j];
                } /* end idx loop */

                lu_solve(A.entries[row_ptr[row]], x0.m[row], swap[row]);
              } /* end row loop */
            } /* end task */
          } /* end tid loop */
        } /* calc */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

      } while(iter < max_iter);
    } /* end single */
  } /* end of parallel region */
} /* gauss_seidel() */
#endif /* USE_TASKLOOP == 1 */

#if USE_TASKLOOP == 0
/*******************************************************************************
*
*******************************************************************************/
static void symm_gauss_seidel(BSMLinSys *linsys, CommMap *parallel_data,
                              Matrix x0, int max_iter)
{
  BlockSparseMatrix A = linsys->A;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  int *row_ptr = A.row_ptr;
  int *col_idx = A.col_index;
  const int brows = A.block_size_row;
  const int bcols = A.block_size_col;
  const int num_rows = A.num_rows;
  const int rhs_rows = rhs.rows;

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

  Matrix help  = generateMatrix(rhs.rows, rhs.cols);

# pragma omp parallel default(none) shared(A, col_idx, row_ptr, rhs, help, \
                                           parallel_data, swap, max_iter, x0)
  {
    int nthreads = omp_get_num_threads();
    int iter = 0;
    int i, j, k;
    CHECK(nthreads > 0);

    /* for loops with one chunk per thread */
    const int chunkmin = num_rows / nthreads;
    const int rest = num_rows % nthreads;

#   pragma omp single
    {
      /*------------------------------------------------------------------------
      | perform a symmetric Gauss-Seidel iteration
      ------------------------------------------------------------------------*/
      do
      {
        iter++;

        /*----------------------------------------------------------------------
        | Forward sweep
        ----------------------------------------------------------------------*/
        /* copy old solution */
#       pragma omp taskgroup
        {
          for(k = 0; k < rhs_rows; k += chunksize)
          {
#           pragma omp task
            {
              for(i = k; i < MIN(k + chunksize, rhs_rows); i++)
              {
                for(j = 0; j < rhs.cols; j++)
                  help.m[i][j] = x0.m[i][j];
              }
            }
          }
        } /* copy */

        /* calculate new solution start...stop-1 */
#       pragma omp taskgroup
        {
          for(int tid = 0; tid < nthreads; tid++)
          {
            int mychunksize = chunkmin + (tid < rest);
            int start = tid * chunkmin + MIN(tid, rest);
            int stop = start + mychunksize;

#           pragma omp task
            {
              for(int row = start; row < stop; row++)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(int idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  double *x0_row = NULL;
                  int col = col_idx[idx];

                  if(col >= start && col < stop)
                    x0_row = x0.m[col]; /* use own new solution */
                  else
                    x0_row = help.m[col]; /* use old solution */

                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= A.entries[idx].m[i][j] * x0_row[j];
                } /* end idx loop */

                lu_solve(A.entries[row_ptr[row]], x0.m[row], swap[row]);
              } /* end row loop */
            } /* end task */
          } /* end tid loop */
        } /* calc */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

        /*----------------------------------------------------------------------
        | Backward sweep
        ----------------------------------------------------------------------*/
        /* copy old solution */
#       pragma omp taskgroup
        {
          for(k = 0; k < rhs_rows; k += chunksize)
          {
#           pragma omp task
            {
              for(i = k; i < MIN(k + chunksize, rhs_rows); i++)
              {
                for(j = 0; j < rhs.cols; j++)
                  help.m[i][j] = x0.m[i][j];
              }
            }
          }
        } /* copy */

        /* calculate new solution stop-1...start */
#       pragma omp taskgroup
        {
          for(int tid = 0; tid < nthreads; tid++)
          {
            int mychunksize = chunkmin + (tid < rest);
            int start = tid * chunkmin + MIN(tid, rest);
            int stop = start + mychunksize;

#           pragma omp task
            {
              for(int row = stop - 1; row >= start; row--)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(int idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  double *x0_row = NULL;
                  int col = col_idx[idx];

                  if(col >= start && col < stop)
                    x0_row = x0.m[col]; /* use own new solution */
                  else
                    x0_row = help.m[col]; /* use old solution */

                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= A.entries[idx].m[i][j] * x0_row[j];
                } /* end idx loop */

                lu_solve(A.entries[row_ptr[row]], x0.m[row], swap[row]);
              } /* end row loop */
            } /* end task */
          } /* end tid loop */
        } /* calc */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

      } while(iter < max_iter);
    } /* end single */
  } /* end of parallel region */
} /* symm_gauss_seidel() */
#endif /* USE_TASKLOOP == 0 */

#if USE_TASKLOOP == 1
/*******************************************************************************
*
*******************************************************************************/
static void symm_gauss_seidel(BSMLinSys *linsys, CommMap *parallel_data,
                              Matrix x0, int max_iter)
{
  BlockSparseMatrix A = linsys->A;
  Matrix rhs = linsys->rhs;
  int **swap = linsys->swap;

  int *row_ptr = A.row_ptr;
  int *col_idx = A.col_index;
  const int brows = A.block_size_row;
  const int bcols = A.block_size_col;
  const int num_rows = A.num_rows;
  const int rhs_rows = rhs.rows;

  CHECK(linsys->decomposition_state == BSM_DIAG_LU_PIVOT);
  CHECK(linsys->diag_first);

  Matrix help  = generateMatrix(rhs.rows, rhs.cols);

# pragma omp parallel default(none) shared(A, col_idx, row_ptr, rhs, help, \
                                           parallel_data, swap, max_iter, x0)
  {
    int nthreads = omp_get_num_threads();
    int iter = 0;
    int i, j;

    CHECK(nthreads > 0);

    /* one chunk per thread */
    const int chunkmin = num_rows / nthreads;
    const int rest = num_rows % nthreads;

#   pragma omp single
    {
      /*------------------------------------------------------------------------
      | perform a symmetric Gauss-Seidel iteration
      ------------------------------------------------------------------------*/
      do
      {
        iter++;

        /*----------------------------------------------------------------------
        | Forward sweep
        ----------------------------------------------------------------------*/
        /* copy old solution */
#       pragma omp taskloop GRAINSIZECLAUSE
        for(int row = 0; row < rhs_rows; row++)
        {
          for(j = 0; j < rhs.cols; j++)
            help.m[row][j] = x0.m[row][j];
        }

        /* calculate new solution start...stop-1 */
#       pragma omp taskgroup
        {
          for(int tid = 0; tid < nthreads; tid++)
          {
            int mychunksize = chunkmin + (tid < rest);
            int start = tid * chunkmin + MIN(tid, rest);
            int stop = start + mychunksize;

#           pragma omp task
            {
              for(int row = start; row < stop; row++)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(int idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  double *x0_row = NULL;
                  int col = col_idx[idx];

                  if(col >= start && col < stop)
                    x0_row = x0.m[col]; /* use own new solution */
                  else
                    x0_row = help.m[col]; /* use old solution */

                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= A.entries[idx].m[i][j] * x0_row[j];
                } /* end idx loop */

                lu_solve(A.entries[row_ptr[row]], x0.m[row], swap[row]);
              } /* end row loop */
            } /* end task */
          } /* end tid loop */
        } /* calc */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait


        /*----------------------------------------------------------------------
        | Backward sweep
        ----------------------------------------------------------------------*/
        /* copy old solution */
#       pragma omp taskloop GRAINSIZECLAUSE
        for(int row = 0; row < rhs_rows; row++)
        {
          for(j = 0; j < rhs.cols; j++)
            help.m[row][j] = x0.m[row][j];
        }

        /* calculate new solution stop-1...start */
#       pragma omp taskgroup
        {
          for(int tid = 0; tid < nthreads; tid++)
          {
            int mychunksize = chunkmin + (tid < rest);
            int start = tid * chunkmin + MIN(tid, rest);
            int stop = start + mychunksize;

#           pragma omp task
            {
              for(int row = stop - 1; row >= start; row--)
              {
                for(i = 0; i < rhs.cols; i++)
                  x0.m[row][i] = rhs.m[row][i];

                for(int idx = row_ptr[row] + 1; idx < row_ptr[row + 1]; idx++)
                {
                  double *x0_row = NULL;
                  int col = col_idx[idx];

                  if(col >= start && col < stop)
                    x0_row = x0.m[col]; /* use own new solution */
                  else
                    x0_row = help.m[col]; /* use old solution */

                  for(i = 0; i < brows; i++)
                    for(j = 0; j < bcols; j++)
                      x0.m[row][i] -= A.entries[idx].m[i][j] * x0_row[j];
                } /* end idx loop */

                lu_solve(A.entries[row_ptr[row]], x0.m[row], swap[row]);
              } /* end row loop */
            } /* end task */
          } /* end tid loop */
        } /* calc */

        /* communicate new solution */
#       pragma omp task
        exchange_matrix(parallel_data, x0);
#       pragma omp taskwait

      } while(iter < max_iter);
    } /* end single */
  } /* end of parallel region */
} /* symm_gauss_seidel() */
#endif /* USE_TASKLOOP == 1 */
