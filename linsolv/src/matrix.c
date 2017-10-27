/*
 * matrix.c
 */

#include "matrix.h"

#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include "util.h"  /* memory allocation */

/*******************************************************************************
*
*******************************************************************************/
Matrix generateMatrix(const int rows, const int cols)
{
  int i;
  Matrix matrix;
  matrix.cols = cols;
  matrix.rows = rows;

  /*----------------------------------------------------------------------------
  | Memory for row pointers
  ----------------------------------------------------------------------------*/
  matrix.m = check_malloc(rows * sizeof(double*));

  /*----------------------------------------------------------------------------
  | Memory for rows
  ----------------------------------------------------------------------------*/
  for(i = 0; i < rows; i++)
    matrix.m[i] = check_malloc(cols * sizeof(double));

  clearMatrix(matrix);

  return matrix;
} /* generateMatrix() */

/*******************************************************************************
*
*******************************************************************************/
extern void deleteMatrix(Matrix matrix)
{
  int i;

  for(i = 0; i < matrix.rows; i++)
    check_free(matrix.m[i]);

  check_free(matrix.m);
} /* deleteMatrix() */

/*******************************************************************************
*
*******************************************************************************/
void clearMatrix(Matrix A)
{
  int i, j;

  for(i = 0; i < A.rows; i++)
    for(j = 0; j < A.cols; j++)
      A.m[i][j] = 0.0;
} /* clearMatrix() */

/*******************************************************************************
*
*******************************************************************************/
BlockSparseMatrix generate_block_sparse_matrix(const int num_rows,
                                               const int num_max_entries,
                                               const int block_size_row,
                                               const int block_size_col)
{
  int i;
  BlockSparseMatrix A;

  A.num_rows       = num_rows;
  A.row_ptr        = check_calloc(num_rows + 1, sizeof(int));
  A.col_index      = check_calloc(num_max_entries, sizeof(int));
  A.entries        = check_malloc(num_max_entries * sizeof(Matrix));
  A.block_size_row = block_size_row;
  A.block_size_col = block_size_col;
  A.size           = num_max_entries;

  for(i = 0; i < num_max_entries; i++)
    A.entries[i]  = generateMatrix(block_size_row, block_size_col);

  return A;
} /* generate_block_sparse_matrix() */

/*******************************************************************************
*
*******************************************************************************/
void delete_block_sparse_matrix(BlockSparseMatrix A)
{
  int i;

  for(i = 0; i < A.size; i++)
    deleteMatrix(A.entries[i]);

  check_free(A.entries);

  check_free(A.row_ptr);

  check_free(A.col_index);
} /* delete_block_sparse_matrix() */

/*******************************************************************************
*
*******************************************************************************/
int check_diag_entry_is_first(BlockSparseMatrix *bsm)
{
  for(int row = 0; row < bsm->num_rows; row++)
    if(bsm->col_index[bsm->row_ptr[row]] != row)
      return 0;

  return 1;
} /* check_diag_entry_is_first() */

/*******************************************************************************
*
*******************************************************************************/
void print_colum_checksums_6(Matrix mat, int nrows)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc = 0;
  int myrank = 0;
  double cksums[6] = {0.0};
  double *recv_cksums = NULL;

  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nproc);

  if(myrank == 0)
  {
    recv_cksums = check_malloc(nproc * 6 * sizeof(double));
  }

  /* calc local checksums */
  for(int i = 0; i < MIN(nrows, mat.rows); i++)
    for(int j = 0; j < MIN(6, mat.cols); j++)
      cksums[j] += mat.m[i][j];

  MPI_Gather(cksums, 6, MPI_DOUBLE, recv_cksums, 6, MPI_DOUBLE, 0, comm);

  /* calc and print global checksums */
  if(myrank == 0)
  {
    double global_cksums[6] = {0.0};

    for(int i = 0; i < nproc; i++)
      for(int j = 0; j < 6; j++)
        global_cksums[j] += recv_cksums[i * 6 + j];

    printf("[0] Column checksums:\n");
    for(int i = 0; i < MIN(6, mat.cols); i++)
      printf("[0]  col %d: %e\n", i, global_cksums[i]);
    fflush(stdout);

    check_free(recv_cksums);
  }
} /* print_colum_checksums_6() */
