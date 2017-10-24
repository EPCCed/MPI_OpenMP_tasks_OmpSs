#ifndef SparseSymmetric

#define SparseSymmetric 1

#include "SparseMatrices.h"

extern void DesymmetrizeSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst);

extern void RestrictUpperSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst);

extern void ProdSymSparseMatrixVector (SparseMatrix spr, double *vec, double *res);

extern void ProdSymSparseMatrixVector2 (SparseMatrix spr, double *vec, double *res);

extern void ProdSymSparseMatrixVector3 (SparseMatrix spr, double *vec, double *res);

#endif
