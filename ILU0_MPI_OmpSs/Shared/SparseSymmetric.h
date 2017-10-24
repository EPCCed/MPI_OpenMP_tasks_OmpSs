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
#ifndef SparseSymmetric

#define SparseSymmetric 1

#include "SparseMatrices.h"

extern void DesymmetrizeSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst);

extern void RestrictUpperSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst);

extern void ProdSymSparseMatrixVector (SparseMatrix spr, double *vec, double *res);

extern void ProdSymSparseMatrixVector2 (SparseMatrix spr, double *vec, double *res);

extern void ProdSymSparseMatrixVector3 (SparseMatrix spr, double *vec, double *res);

#endif
