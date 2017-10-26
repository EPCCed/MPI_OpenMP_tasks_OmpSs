#ifndef SparseHarwellBoeing

#define SparseHarwellBoeing 1

#include "SparseMatrices.h"

extern void ReadSparseIndexValuesHB (char *nameFile, int *vpos, double *vval, int FtoC);

extern void ReadSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC);

extern void CreateSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC);

extern void CreateSparseMatrixHB2 (char *nameFile, ptr_SparseMatrix spr, int FtoC);

extern void WriteSparseMatrixHB (char *nameFile, SparseMatrix spr, int numI, int numR, 
													char *tag, int CtoF);

#endif
