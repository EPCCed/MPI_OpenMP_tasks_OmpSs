/*
 * matrix.h
 */

#ifndef SRC_MATRIX_H_
#define SRC_MATRIX_H_

/*******************************************************************************
* Basic definition of matrix structures and functions to allocate memory
* and free memory. Implemented structures are general dense matrices given
* by the data structure Matrix. A data structure for block sparse matrices
* is given by BlockSparseMatrix. These matrices are given by the Block
* Compressed Sparse Row (BCSR) format.
*******************************************************************************/

/*******************************************************************************
* General definition of a dense matrix with double entries.
*******************************************************************************/
typedef struct
{
  double** m; /* Array for matrix entries of size rows * cols */
  int cols;   /* Number of columns                            */
  int rows;   /* Number of rows                               */
} Matrix;

/*******************************************************************************
* Structure defining a block sparse matrix in Block Compressed Sparse Row
* (BCSR) format.
*******************************************************************************/
typedef struct
{
  int block_size_row;   /* Number of rows of each block              */
  int block_size_col;   /* Number of columns of each block           */

  int num_rows;         /* Number of rows of the block sparse matrix */
  int *row_ptr;         /* CSR format: Index of first entry in each row */
  int *col_index;       /* CSR format: Corresponding column index    */

  int size;             /* Size of block sparse matrix               */

  Matrix *entries;      /* Pointer to entries in block sparse matrix */
} BlockSparseMatrix;

/*******************************************************************************
* Generates Matrix. Memory for Matrix.m (rows * cols doubles) will be allocated
* and initialized to zero.
*
* @param[in] rows Number of rows of Matrix.
* @param[in] cols Number of columns of Matrix.
* @return Matrix with dimension rows times cols of type double.
*******************************************************************************/
Matrix generateMatrix(const int rows, const int cols);

/*******************************************************************************
* Frees memory allocated for a Matrix, i.e.\ the memory allocated for the double
* array matrix.m is freed.
*
* @param[in] matrix Matrix for which allocated memory is freed
*******************************************************************************/
void deleteMatrix(Matrix matrix);

/*******************************************************************************
* Clears the matrix A, that is all entries are set to zero.
*
* @param[in,out] A Matrix which is cleared to zero.
*******************************************************************************/
void clearMatrix(Matrix A);

/*******************************************************************************
* Allocates memory for a block sparse matrix in Block Compressed Sparse Row
* (BCSR) format for a given number of rows. Each entry in the matrix is a
* Matrix with given a fixed number of rows and columns.
*
* @param[in] num_rows Number of rows of the block sparse matrix
* @param[in] num_max_entries Maximum number of entries, for which memory is
*            allocated
* @param[in] block_size_row Number of row of each entry of the sparse matrix
* @param[in] block_size_col Number of columns of each entry of the sparse matrix
* @return A clear block sparse matrix for saving a matrix with the defined
*         number of rows and maximum entries
*******************************************************************************/
BlockSparseMatrix generate_block_sparse_matrix(const int num_rows,
                                               const int num_max_entries,
                                               const int block_size_row,
                                               const int block_size_col);

/*******************************************************************************
* Frees memory which was allocated for a block sparse matrix in Block Compressed
* Sparse Row (BCSR) format.
*
* @param[in] A Block sparse matrix for which memory is freed.
*******************************************************************************/
void delete_block_sparse_matrix(BlockSparseMatrix A);

/*******************************************************************************
* Test if diagonal entry is stored first for each column or not.
*******************************************************************************/
int check_diag_entry_is_first(BlockSparseMatrix *bsm);

/*******************************************************************************
* Calculate and print simple checksum for matrix colums (up to 6 colums and
* nrows rows).
*******************************************************************************/
void print_colum_checksums_6(Matrix mat, int nrows);

#endif /* SRC_MATRIX_H_ */
