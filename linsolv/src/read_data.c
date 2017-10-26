/*
 * read_data.c
 */

#include "read_data.h"

#include <stddef.h> /* NULL */
#include <stdlib.h> /* EXIT_FAILURE */
#include <stdio.h>
#include <netcdf.h>
#include "util.h"
#include "matrix.h"
#include "setup_comm.h"

#define check_ncerr(status, name) \
  if((status != NC_NOERR)) \
  { \
    printf("Error: %s (%s) [%s:%i]\n", nc_strerror(status), name, \
           __FILE__, __LINE__); \
    exit(EXIT_FAILURE); \
  }

static int open_ncfile(char *filename);
static void close_ncfile(int ncid);
static int read_ncdimension(int ncid, const char *name);
static void read_ncvariable(int ncid, const char *name, int type, void *data);

/*******************************************************************************
*
*******************************************************************************/
static int open_ncfile(char *filename)
{
  int ncid;

  check_ncerr(nc_open(filename, NC_NOWRITE, &ncid), filename);

  return ncid;
} /* open_ncfile() */

/*******************************************************************************
*
*******************************************************************************/
static void close_ncfile(int ncid)
{
  check_ncerr(nc_close(ncid), "-");
} /* close_ncfile() */

/*******************************************************************************
*
*******************************************************************************/
static int read_ncdimension(int ncid, const char *name)
{
  int dimid;
  size_t dimlen;

  check_ncerr(nc_inq_dimid(ncid, name, &dimid), name);
  check_ncerr(nc_inq_dimlen(ncid, dimid, &dimlen), name);

  return dimlen;
} /* close_ncfile() */

/*******************************************************************************
*
*******************************************************************************/
static void read_ncvariable(int ncid, const char *name, int type, void *data)
{
  int varid;

  check_ncerr(nc_inq_varid(ncid, name, &varid), name);

  switch(type)
  {
    case NC_INT:
      check_ncerr(nc_get_var_int(ncid, varid, (int *)data), name);
      break;
    case NC_DOUBLE:
      check_ncerr(nc_get_var_double(ncid, varid, (double *)data), name);
      break;
    default:
      printf("Error: unknown type (%s)\n", STR(type));
      exit(EXIT_FAILURE);
  }
} /* close_ncfile() */

/*******************************************************************************
*
*******************************************************************************/
void read_comap(char *filename, CommMap *comap)
{
  /* open the file */
  int ncid = open_ncfile(filename);
  DBG_MSG("  reading file %s\n", filename);

  /* read dimensions */
  int ndom = read_ncdimension(ncid, "ndomains");
  int ncommdom = read_ncdimension(ncid, "ncommdomains");
  int nown = read_ncdimension(ncid, "nownpoints");
  int nadd = read_ncdimension(ncid, "naddpoints");

  /* check */
  CHECK(ndom == comap->ndomains);

  /* alloc memory */
  int *addpoint_owner = check_calloc(nadd, sizeof(int));
  int *addpoint_id = check_calloc(nadd, sizeof(int));
  int *commpartner = check_calloc(ncommdom, sizeof(int));
  int *sendcount = check_calloc(ndom, sizeof(int));
  int *recvcount = check_calloc(ndom, sizeof(int));

  /* read data */
  read_ncvariable(ncid, "addpoint_owner", NC_INT, addpoint_owner);
  read_ncvariable(ncid, "addpoint_idx", NC_INT, addpoint_id);
  read_ncvariable(ncid, "commpartner", NC_INT, commpartner);
  read_ncvariable(ncid, "sendcount", NC_INT, sendcount);
  read_ncvariable(ncid, "recvcount", NC_INT, recvcount);

  /* copy read data to comap structure */
  CHECK(comap != NULL);
  comap->ncommdomains = ncommdom;
  comap->nownpoints = nown;
  comap->naddpoints = nadd;
  comap->addpoint_owner = addpoint_owner;
  comap->addpoint_id = addpoint_id;
  comap->commpartner = commpartner;
  comap->sendcount = sendcount;
  comap->recvcount = recvcount;

  /* close file. */
  close_ncfile(ncid);
} /* read_comap() */

/*******************************************************************************
*
*******************************************************************************/
void read_block_sparse_matrix(char *filename, BlockSparseMatrix *bsmat)
{
  int num_rows = 0;
  int bsmatrix_size = 0;
  int block_size_row = 0;
  int block_size_col = 0;
  int dim_bsmatrix_1d = 0;

  int *row_ptr = NULL;
  int *col_index = NULL;
  Matrix *entries = NULL;
  double *bsmatrix_1d = NULL;

  /* open the file */
  int ncid = open_ncfile(filename);
  DBG_MSG("  reading file %s\n", filename);

  /* read dimensions */
  num_rows = read_ncdimension(ncid, "num_rows");
  bsmatrix_size = read_ncdimension(ncid, "bsmatrix_size");
  block_size_row = read_ncdimension(ncid, "block_size_row");
  block_size_col = read_ncdimension(ncid, "block_size_col");
  dim_bsmatrix_1d = read_ncdimension(ncid, "dim_bsmatrix_1d");

  DBG_MSG("  dimensions: num_rows: %d, size: %d, "
          "block_size_row: %d, block_size_col: %d\n",
          num_rows, bsmatrix_size, block_size_row, block_size_col);

  /* check dimensions */
  CHECK(num_rows > 0);
  CHECK(bsmatrix_size > 0);
  CHECK(block_size_row > 0);
  CHECK(block_size_col > 0);
  CHECK(dim_bsmatrix_1d > 0);

  /* alloc memory */
  row_ptr = check_calloc(num_rows + 1, sizeof(int));
  col_index = check_calloc(bsmatrix_size, sizeof(int));
  entries = check_malloc(bsmatrix_size * sizeof(Matrix));
  for(int i = 0; i < bsmatrix_size; i++) /* alloc submatrices */
    entries[i] = generateMatrix(block_size_row, block_size_col);
  bsmatrix_1d = check_malloc(dim_bsmatrix_1d * sizeof(double));

  /* read data */
  read_ncvariable(ncid, "col_index", NC_INT, col_index);
  read_ncvariable(ncid, "num_entries_row", NC_INT, row_ptr);
  read_ncvariable(ncid, "bsmatrix_1d", NC_DOUBLE, bsmatrix_1d);

  /* copy data from 1d matrix to sub matrices */
  int idx = 0;
  for(int i = 0; i < bsmatrix_size; i++)
    for(int row = 0; row < block_size_row; row++)
      for(int col = 0; col < block_size_col; col++)
        entries[i].m[row][col] = bsmatrix_1d[idx++];

  check_free(bsmatrix_1d);

  /* copy read data to bsmatrix structure */
  CHECK(bsmat != NULL);
  bsmat->num_rows = num_rows;
  bsmat->size = bsmatrix_size;
  bsmat->block_size_row = block_size_row;
  bsmat->block_size_col = block_size_col;
  bsmat->row_ptr = row_ptr;
  bsmat->col_index = col_index;
  bsmat->entries = entries;

  /* close file. */
  close_ncfile(ncid);
} /* read_block_sparse_matrix() */


/*******************************************************************************
*
*******************************************************************************/
void read_matrix(char *filename, Matrix *mat)
{
  int rows = 0;
  int cols = 0;
  int dim_matrix_1d = 0;

  double **matrix = NULL;
  double *matrix_1d = NULL;

  /* open the file */
  int ncid = open_ncfile(filename);
  DBG_MSG("  reading file %s\n", filename);

  /* read dimensions */
  rows = read_ncdimension(ncid, "size_row");
  cols = read_ncdimension(ncid, "size_col");
  dim_matrix_1d = read_ncdimension(ncid, "dim_matrix_1d");
  DBG_MSG("  dimensions: rows: %d, cols: %d\n", rows, cols);

  /* check dimensions */
  CHECK(rows > 0);
  CHECK(cols > 0);
  CHECK(dim_matrix_1d > 0);
  CHECK(dim_matrix_1d == rows * cols);

  /* alloc memory */
  matrix = check_malloc(rows * sizeof(double*));
  for(int i = 0; i < rows; i++)
    matrix[i] = check_malloc(cols * sizeof(double));
  matrix_1d = check_malloc(dim_matrix_1d * sizeof(double));

  /* read data */
  read_ncvariable(ncid, "matrix_1d", NC_DOUBLE, matrix_1d);

  /* copy data from 1d matrix to matrix */
  int idx = 0;
  for(int row = 0; row < rows; row++)
    for(int col = 0; col < cols; col++)
      matrix[row][col] = matrix_1d[idx++];

  check_free(matrix_1d);

  /* copy read data to matrix structure */
  CHECK(mat != NULL);
  mat->rows = rows;
  mat->cols = cols;
  mat->m = matrix;

  /* close file. */
  close_ncfile(ncid);
} /* read_matrix() */
