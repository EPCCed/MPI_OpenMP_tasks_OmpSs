#ifndef SparseMatrixTip

#define SparseMatrixTip 1

#include "SparseVectors.h"

typedef struct
	{
		int dim1, dim2;
		int *vptr;
		int *vpos;
		double *vval;
	} SparseMatrix, *ptr_SparseMatrix;

extern void CreateSparseMatrix (ptr_SparseMatrix spr, int numR, int numC, 
															 int numE, int msr);

extern void ReallocSparseMatrix (ptr_SparseMatrix spr, int msr);

extern void ConvertSparseMatrix (SparseMatrix spr, int lngs);

extern void CopySparseMatrices (SparseMatrix src, ptr_SparseMatrix dst, int lngs);

extern void PrintSparseMatrix (SparseMatrix spr, int CorF);

extern void FPrintSparseMatrix (FILE *f, SparseMatrix spr, int CorF);

extern void RemoveSparseMatrix (ptr_SparseMatrix spr);

extern void GetDiagonalSparseMatrix (SparseMatrix spr, double *diag);

extern void GetDiagonalSparseMatrixDspls (SparseMatrix spr, double *diag, int dspls);

extern void HeapSortSparseMatrix (SparseMatrix spr, SparseVector_Sort_func funcSort);

extern void PermuteColsSparseMatrix (SparseMatrix spr, int *perm);

extern void PermuteColsWithNegSparseMatrix (SparseMatrix spr, int *perm);

extern void GetGraphWithDiagSparseMatrix (SparseMatrix spr, int *vptr, int *vpos);

extern void GetGraphSparseMatrix (SparseMatrix spr, int *vptr, int *vpos);

extern void GetGraphSparseMatrix2 (SparseMatrix spr, int *vptr, int *vpos);

extern void GetGraphSparseMatrix3 (SparseMatrix spr, int indexS, int *vptr, int *vpos, int indexG);

extern void ProdSparseMatrixVector (SparseMatrix spr, double *vec, double *res);

extern void ProdSparseMatrixVector2 (SparseMatrix spr, double *vec, double *res);

extern void ProdSparseMatrixVector3 (SparseMatrix spr, double *vec, double *res);


extern void AddSparseMatrices (SparseMatrix src1, int *perm1, int *iperm1,
												SparseMatrix src2, int *perm2, int *iperm2,
                        ptr_SparseMatrix dst);

#endif

