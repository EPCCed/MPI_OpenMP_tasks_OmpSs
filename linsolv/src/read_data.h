/*
 * read_data.h
 */

#ifndef SRC_READ_DATA_H_
#define SRC_READ_DATA_H_

#include "matrix.h"
#include "setup_comm.h"

/*******************************************************************************
*
*******************************************************************************/
void read_comap(char *filename, CommMap *comap);

/*******************************************************************************
* Read block sparse matrix from NetCDF file.
*******************************************************************************/
void read_block_sparse_matrix(char *filename, BlockSparseMatrix *bsmat);

/*******************************************************************************
* Read matrix from NetCDF file.
*******************************************************************************/
void read_matrix(char *filename, Matrix *mat);


#endif /* SRC_READ_DATA_H_ */
