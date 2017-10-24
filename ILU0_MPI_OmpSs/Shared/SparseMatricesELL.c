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
#include <stdio.h>
#include <stdlib.h>
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectors.h"
#include "SparseMatrices.h"
#include "SparseMatricesELL.h"

void CreateSparseMatrixELL (ptr_SparseMatrixELL ell, int numR, int numC, int numEpR)
	{
		ell->dim1 = numR; ell->dim2 = numC; ell->numEpR = numEpR;
		CreateMatrixInts    (&ell->mpos, numR, numEpR);
		CreateMatrixDoubles (&ell->mval, numR, numEpR);
	}

void RemoveSparseMatrixELL (ptr_SparseMatrixELL ell)
	{
		RemoveMatrixDoubles (&ell->mval);
		RemoveMatrixInts (&ell->mpos);
		ell->dim1 = -1; ell->dim2 = -1; ell->numEpR = -1;
	}

void PrintSparseMatrixELL (SparseMatrixELL ell, int CorF)
	{
		int i, j;

		printf ("Values: \n");
		for (i=0; i<ell.dim1; i++)
			{ printf (" Row %d --> ", i+CorF);
			for (j=0; j<ell.numEpR; j++)
				printf ("(%d,%f) ", ell.mpos[i][j], ell.mval[i][j]); 
			printf ("\n");	}
		printf ("\n");
	}

void ConvertSparseMatrixCSR2ELL (SparseMatrix spr, ptr_SparseMatrixELL ell) 
	{
		int i, max;
		int *vlen = NULL;

		CreateInts (&vlen, spr.dim1);
		ComputeLengthfromHeader (spr.vptr, vlen, spr.dim1);

		max = vlen[0];
		for (i=1; i<spr.dim1; i++)
			if (max < vlen[i]) max = vlen[i]; 

		CreateSparseMatrixELL (ell, spr.dim1, spr.dim2, max);
		InitInts (ell->mpos[0], spr.dim1*max, 0, 0);
		InitDoubles (ell->mval[0], spr.dim1*max, 0.0, 0.0);
		for (i=0; i<spr.dim1; i++)
			{
				CopyInts (spr.vpos+spr.vptr[i], ell->mpos[i], vlen[i]);
				CopyDoubles (spr.vval+spr.vptr[i], ell->mval[i], vlen[i]);
			}

		RemoveInts (&vlen);
	}

void ProdSparseMatrixELLVector (SparseMatrixELL ell, double *vec, double *res)
	{
		int i, j;
		double aux;

		for (i=0; i<ell.dim1; i++)
			{
				aux = 0.0;
				for (j=0; j<ell.numEpR; j++)
					if (ell.mval[i][j] != 0.0)
						aux += ell.mval[i][j] * vec[ell.mpos[i][j]];
				res[i] += aux;
			}
	}

