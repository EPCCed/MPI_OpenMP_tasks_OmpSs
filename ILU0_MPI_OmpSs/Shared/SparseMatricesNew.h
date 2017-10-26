/*
 * Copyright (C) 2016 UJI
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License in COPYING.LGPL for more details.
 */
#ifndef SparseMatrixNewTip

#define SparseMatrixNewTip 1

#include "SparseVectorsNew.h"

typedef struct
	{
		int dim1, dim2;
		int *vptr;
		int *vpos;
		double *vval;
	} SparseMatrix, *ptr_SparseMatrix;

/*******************************************************************/

// This routine creates a sparseMatrix from the next parameters
// * numR defines the number of rows
// * numC defines the number of columns
// * numE defines the number of nonzero elements
// * msr indicates if the MSR is the format used to the sparse matrix
// If msr is actived, numE doesn't include the diagonal elements
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void CreateSparseMatrix (ptr_SparseMatrix p_spr, int index, int numR, int numC, int numE, int msr);

// This routine adjusts the memory size of a sparseMatrix 
extern void ReallocSparseMatrix (ptr_SparseMatrix p_spr);

// This routine liberates the memory related to matrix spr
extern void RemoveSparseMatrix (ptr_SparseMatrix spr);

/*******************************************************************/

// This routine prints on the screen the contents of spr in a simple way.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void PrintSparseMatrix (SparseMatrix spr, int index);

// This routine writes to the file f the contents of spr in a simple way.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void FPrintSparseMatrix (FILE *f, SparseMatrix spr, int index);

// This routine prints on the screen the contents of spr, such that the data
// could be read by MATLAB.	The parameter index indicates if 0-indexing 
// or 1-indexing is used. The values fI1 and fI2 refer to the way in which
// the indices are printed, while the fD1 and fD2 are related to the values.
extern void PrintMatlabFSparseMatrix (SparseMatrix spr, int index, int fI1, int fI2, int fD1, int fD2);

// This routine writes on the file f the contents of spr, such that the data
// could be read by MATLAB.	The parameter index indicates if 0-indexing 
// or 1-indexing is used. The values fI1 and fI2 refer to the way in which
// the indices are printed, while the fD1 and fD2 are related to the values.
extern void FPrintMatlabFSparseMatrix (FILE *f, SparseMatrix spr, int index, 
																				int fI1, int fI2, int fD1, int fD2);

// This routine writes on the file filename the contents of spr, such that the data
// could be read by MATLAB.	The parameter index indicates if 0-indexing 
// or 1-indexing is used. The values fI1 and fI2 refer to the way in which
// the indices are printed, while the fD1 and fD2 are related to the values.
extern void WriteMatlabFSparseMatrix (char *filename, SparseMatrix spr, int index, 
																				int fI1, int fI2, int fD1, int fD2);

// This routine writes different files to store the contents of spr.
// The name of these files begins with the string prefix.
// The parameter index indicates if 0-indexing or 1-indexing is used. 
// The values fI1 and fI2 refer to the way in which the indices are printed, 
// while the fD1 and fD2 are related to the values.
extern void WriteMatlabFSparseMatrix2 (char *prefix, SparseMatrix spr, int index, 
																				int fI1, int fI2, int fD1, int fD2);

/*******************************************************************/

// If it is necessary, this routine changes the indexing of spr,
// from 0-indexing to 1-indexing or in the reverse case.
// The parameter lngs defines the change, where its first bit
// defines the indexing at the input and its second bit defines
// the deseared indexing at the output, thus
// * if is equal to 1, changes to 1 to 0 indexing.
// * if is equal to 2, changes to 0 to 1 indexing.
extern void ConvertSparseMatrix (SparseMatrix spr, int lngs);

/*******************************************************************/

// This routine sorts the rows of spr, according to the criteria defined in funcSort.
// This operation is made row by row.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void HeapSortSparseMatrix (SparseMatrix spr, int index, SparseVector_Sort_func funcSort);

// This routine permuts the columns of spr, if perm is not null,
// or it only sort the rows according to the criteria defined in funcSort.
// The parameter index indicates if 0-indexing or 1-indexing is used,
// for both, the matrix spy and the vector perm.
extern void PermuteColsSparseMatrix (SparseMatrix spr, int index, int *perm);

// This routine permuts the columns of spr, if perm is not null,
// or it only sort the rows according to the criteria defined in funcSort.
// The perm contains negative indices to mark the columns to be removed.
// The parameter index indicates if 0-indexing or 1-indexing is used,
// for both, the matrix spy and the vector perm.
extern void PermuteColsWithNegSparseMatrix (SparseMatrix spr, int index, int *perm);

/*******************************************************************/

// This routine extracts the diagonal of spr on the vector diag.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void GetDiagonalSparseMatrix (SparseMatrix spr, int index, double *diag);

// This routine extracts the diagonal of spr on the vector diag.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void GetDiagonalSparseMatrixDspls (SparseMatrix spr, int index, double *diag, int dspls);

// This routine computes the graph related to the symmetric sparse matrix spr, in which 
// the diagonals are not considered, and the result is stored on the vectors vptr and vpos.
// The parameters indexS and indexG indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrix and the graph.
extern void GetGraphSparseMatrix (SparseMatrix spr, int indexS, int *vptr, int *vpos, int indexG);

// This routine computes the graph related to the symmetric sparse matrix spr, in which 
// the diagonals are considered, and the result is stored on the vectors vptr and vpos.
// The parameters indexS and indexG indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrix and the graph.
extern void GetGraphWithDiagSparseMatrix (SparseMatrix spr, int indexS, int *vptr, int *vpos, int indexG);

/*******************************************************************/

// This routine computes the diagonal scaling { A = diag(rowscal)*A*diag(colscal) }, 
// using sparse vectors.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ScalSparseMatrix (SparseMatrix spr, int index, double *rowscal, double *colscal);

// This routine computes the product { res += spr * vec }, using sparse vectors.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVector (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVector2 (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res = spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVector3 (SparseMatrix spr, int index, double *vec, double *res);

/*******************************************************************/
 
// This routine computes the addition { dst = spr1(perm1,perm1) + spr2(perm2,perm2) }
// The vector iperm1 contains the inverse permutation of perm1, in which the value
// -1 indicates that the correspondient position doesn't exist in perm1.
// The same for iperm2 and perm2.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing 
// or 1-indexing is used to store the different sparse matrices.
// permCols determinates if the vpos vectors of the matrices have to be permuted.
extern void AddSparseMatrices (SparseMatrix src1, int index1, int *perm1, int *iperm1,
																SparseMatrix src2, int index2, int *perm2, int *iperm2,
																ptr_SparseMatrix dst, int index3, int perm_cols);

// This routine computes the number of elements which appears in both sparse matrices.
// The vector iperm1 contains the inverse permutation of perm1, in which the value
// -1 indicates that the correspondient position doesn't exist in perm1.
// The same for iperm2 and perm2.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing 
// or 1-indexing is used to store the different sparse matrices.
extern long int IntersectSparseMatrices (SparseMatrix src1, int index1, int *perm1, int *iperm1,
																					SparseMatrix src2, int index2, int *perm2, int *iperm2,
																					ptr_SparseMatrix dst, int index3);

/*******************************************************************/

#endif

