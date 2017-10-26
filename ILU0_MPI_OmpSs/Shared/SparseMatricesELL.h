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
#ifndef SparseMatrixELLTip

#define SparseMatrixELLTip 1

#include "SparseVectors.h"
#include "SparseMatrices.h"

typedef struct
	{
		int dim1, dim2, numEpR;
		matInts mpos;
		matDoubles mval;
	} SparseMatrixELL, *ptr_SparseMatrixELL;

extern void CreateSparseMatrixELL (ptr_SparseMatrixELL ell, int numR, int numC, int numEpR);

extern void RemoveSparseMatrixELL (ptr_SparseMatrixELL ell);

extern void PrintSparseMatrixELL (SparseMatrixELL ell, int CorF);

extern void ConvertSparseMatrixCSR2ELL (SparseMatrix spr, ptr_SparseMatrixELL ell);

extern void ProdSparseMatrixELLVector (SparseMatrixELL ell, double *vec, double *res);

#endif

