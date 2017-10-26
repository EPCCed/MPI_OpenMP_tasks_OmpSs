/*
 * linsolv.h
 */

#ifndef SRC_LINSOLV_H_
#define SRC_LINSOLV_H_

#include "setup_comm.h"
#include "linsys.h"
#include "matrix.h"

/*******************************************************************************
* Computes LU decompositions of the diagonal entries of the block sparse
* matrix bsmls->A.
*
* @param[in,out] bsmls Given linear system. On output, the diagonal blocks of
*                      bsmls->A are overwritten by LU decompositions with
*                      pivoting and bsmls->swap contains permutation information
*                      of the LU decompositions.
*******************************************************************************/
void linsys_lu_decomp_diag(BSMLinSys *bsmls);

/*******************************************************************************
* This is a wrapper function calling one of the following functions:
* point_implicit(), jacobi(), gauss_seidel(), symm_gauss_seidel().
* The choice of the functions above is done via the variable linear_solver.
* Valid values corresponding to the functions above are
* BSLS_POINT_IMPLICIT, BSLS_JACOBI, BSLS_GAUSS_SEIDEL, BSLS_SYMM_GAUSS_SEIDEL
*
* @param[in] linear_solver Choice of the iterative solution method.
* @param[in] num_sweeps Number of sweeps performed in the iterative solution
*                       method.
* @param[in] linsys     Given linear system. LU decompositions of the diagonal
*                       blocks of linsys->A must have been computed beforehand.
* @param[in,out] approximate_solution On input initial guess for the iterative
*                       solver, on output approximate solution obtained after
*                       given number of sweeps.
* @param[in] comap      Communication map.
*******************************************************************************/
void linsolv(const BSLinSysMethod linear_solver, const int num_sweeps,
             BSMLinSys *linsys, Matrix approximate_solution, CommMap *comap);

#endif /* SRC_LINSOLV_H_ */
