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
#ifndef DoubleVector

#define DoubleVector 1

#include <stdio.h>

extern void CreateDoubles (double **vdbl, int dim);

extern void ReallocDoubles (double **vdbl, int dim);

extern void InitDoubles (double *vdbl, int dim, double frst, double incr);

extern void InitRandDoubles (double *vdbl, int dim, double frst, double last);

extern void CopyDoubles (double *src, double *dst, int dim);

extern void ScaleDoubles (double *vdbl, double scal, int dim);

extern void AxpyDoubles (double alfa, double *vdbl1, double *vdbl2, int dim);

extern void XpayDoubles (double *vdbl1, double alfa, double *vdbl2, int dim);

extern double DotDoubles (double *vdbl1, double *vdbl2, int dim);

extern void GetDoubleFromString (char *string, double *pdbl, int dimC);

extern void GetDoublesFromString (char *string, double *vdbl, int dimN, int dimC);

typedef int (*Doubles_SortPerm_func) (double , double , int, int);

extern int HeapSortPermDoubles_IncrPos (double n1, double n2, int p1, int p2);

extern int HeapSortPermDoubles_DecrPos (double n1, double n2, int p1, int p2);

extern int HeapSortPermDoubles_IncrVal (double n1, double n2, int p1, int p2);

extern int HeapSortPermDoubles_DecrVal (double n1, double n2, int p1, int p2);

extern void HeapSortPermDoubles_Swap (double *n1, double *n2, int *p1, int *p2);

extern void HeapSortPermDoubles_Adjust (double *vdbl, int *prm, int inic, int final, Doubles_SortPerm_func funcSort);

extern void HeapSortPermDoubles (double *vdbl, int *prm, int dim, Doubles_SortPerm_func funcSort);

typedef int (*Doubles_Sort_func) (double , double );

extern int HeapSortDoubles_IncrPos (double n1, double n2);

extern int HeapSortDoubles_DecrPos (double n1, double n2);

extern void HeapSortDoubles_Swap (double *n1, double *n2);

extern void HeapSortDoubles_Adjust (double *vdbl, int inic, int final, Doubles_Sort_func funcSort);

extern void HeapSortDoubles (double *vdbl, int dim, Doubles_Sort_func funcSort);

extern int FindSortDouble (int num, double *vdbl, int dim, Doubles_Sort_func funcSort);

extern void VdivDoubles (double alfa, double *src, double *dst, int dim);

extern void VvecDoubles (double alfa, double *src1, double *src2, double beta, double *dst, int dim);

extern double AddDoubles (double *vdbl, int dim);

extern double AddPermuteDoubles (double *vdbl, int *perm, int index, int dim);

extern double MaxDoubles (double *vdbl, int dim);

extern void CopyPermuteDoubles (double *src, int *perm, int index, double *dst, int dim);

extern void CopyInvPermuteDoubles (double *src, double *dst, int *perm, int index, int dim);

extern void PrintDoubles (double *vdbl, int dim);

extern void FPrintDoubles (FILE *f, double *vdbl, int dim);

extern void PrintFDoubles (double *vdbl, int dim, int f1, int f2);

extern void FPrintFDoubles (FILE *f, double *vdbl, int dim, int f1, int f2);

extern int ReadDoubles (char *file, double **vdbl);

extern void WriteDoubles (char *filename, double *vdbl, int dim);

extern void WriteFDoubles (char *filename, double *vdbl, int dim, int f1, int f2);

extern void RemoveDoubles (double **vdbl);

/****************************************************************************/

typedef double ** matDoubles;

extern void CreateVectorPtrDoubles (matDoubles *mdbl, int numR);

extern void CopyVectorPtrDoubles (matDoubles msrc, matDoubles mdst, int dim);

extern void DesplVectorPtrDoubles (matDoubles mdbl, int dim, int dspl);

extern void RemoveVectorPtrDoubles (matDoubles *mdbl);

extern void CreateMatrixDoubles (matDoubles *mdbl, int numR, int numC);

extern void BuildMatrixDoubles (matDoubles mdbl, double *vec, int numR, int numC);

extern void RemoveMatrixDoubles (matDoubles *mdbl);

#endif

