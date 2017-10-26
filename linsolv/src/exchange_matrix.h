/*
 * exchange_matrix.h
 */

#ifndef SRC_EXCHANGE_MATRIX_H_
#define SRC_EXCHANGE_MATRIX_H_

#include "setup_comm.h"
#include "matrix.h"

/*******************************************************************************
*
*******************************************************************************/
void matrix_comm_init(CommMap *comap, int ncols);

/*******************************************************************************
*
*******************************************************************************/
void matrix_comm_end(void);

/*******************************************************************************
* Communicates the overlapping points in the given matrix between the
* processes defined in the data structure required for parallelization.
*
* @param[in]     comap   Information required for parallelization
* @param[in,out] matrix  Data which is communicated between the different
*                        processes. The information in the additional points
*                        because of parallelization is overwritten by the
*                        information of other processes.
*******************************************************************************/
void exchange_matrix(const CommMap *comap, Matrix matrix);

#endif /* SRC_EXCHANGE_MATRIX_H_ */
