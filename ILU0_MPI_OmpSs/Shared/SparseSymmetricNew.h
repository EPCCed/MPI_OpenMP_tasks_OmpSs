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
#ifndef SparseSymmetricNew

#define SparseSymmetricNew 1

#include "SparseMatricesNew.h"


/*******************************************************************/

// This routine creates de sparse matrix dst from the symmetric matrix spr.
// The parameters indexS and indexG indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices.
extern void DesymmetrizeSparseMatrices (SparseMatrix src, int indexS, ptr_SparseMatrix dst, int indexD);

// This routine creates de sparse matrix dst from the symmetric matrix spr, 
// on which only the upper triangle is stored
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices.
extern void RestrictUpperSparseMatrices (SparseMatrix src, int indexS, ptr_SparseMatrix dst, int indexD);
/*
extern void RestrictUpperSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst);
*/

/*******************************************************************/

// This routine computes the product { res += spr * vec }, using sparse vectors.
// The matrix spr is symmetric, and only one of the replicated value is stored.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSymSparseMatrixVector (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res += spr * vec }.
// The matrix spr is symmetric, and only one of the replicated value is stored.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSymSparseMatrixVector2 (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res = spr * vec }.
// The matrix spr is symmetric, and only one of the replicated value is stored.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSymSparseMatrixVector3 (SparseMatrix spr, int index, double *vec, double *res);

/*******************************************************************/

#endif

