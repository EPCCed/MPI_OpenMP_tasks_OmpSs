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

