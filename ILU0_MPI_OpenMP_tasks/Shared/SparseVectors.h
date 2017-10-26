#ifndef SparseVectorTip

#define SparseVectorTip 1

#include <stdio.h>

typedef struct SparseVector {
	int dim;
	int *vpos;
	double *vval;
	} SparseVector, *ptr_SparseVector;

extern void CreateSparseVector (ptr_SparseVector b, int nz);

extern void CopySparseVectors (SparseVector src, ptr_SparseVector dst);

extern void PrintSparseVectors (SparseVector b);

extern void FPrintSparseVectors (FILE *f, SparseVector b);

extern void RemoveSparseVector (ptr_SparseVector b);

typedef int (*SparseVector_Sort_func) (int , int , double, double);

extern int HeapSortSparseVector_IncrPosWithNeg (int n1, int n2, double x1, 
																								double x2);

extern int HeapSortSparseVector_DecrPosWithNeg (int n1, int n2, double x1, 
																								double x2);

extern int HeapSortSparseVector_IncrPos (int n1, int n2, double x1, double x2);

extern int HeapSortSparseVector_DecrPos (int n1, int n2, double x1, double x2);

extern int HeapSortSparseVector_IncrVal (int n1, int n2, double x1, double x2);

extern int HeapSortSparseVector_DecrVal (int n1, int n2, double x1, double x2);

extern int HeapSortSparseVector_IncrAbsVal (int n1, int n2, double x1, double x2);

extern int HeapSortSparseVector_DecrAbsVal (int n1, int n2, double x1, double x2);

extern void HeapSortSparseVector_Swap (int *n1, int *n2, double *x1, 
																			 double *x2);

extern void HeapSortSparseVector_Adjust (int *vpos, double *vval, int inic, 
																				 int final,
																				 SparseVector_Sort_func funcSort);

extern void HeapSortSparseVector (SparseVector b, 
																	SparseVector_Sort_func funcSort);

extern void GatherSparseVector (double *v, int n, ptr_SparseVector b);

extern void GatherZeroSparseVector (double *v, int n, ptr_SparseVector b);

extern void ScatterSparseVector (SparseVector b, double *v);

extern double DotSparseVector (SparseVector b, double *v);

extern double DotSparseVectors (SparseVector b1, SparseVector b2);

extern void AxpySparseVector (SparseVector b, double alfa, double *v);

extern void ScalSparseVector (SparseVector b, double alfa);

extern void VvecSparseVector (SparseVector b, double *v1, double alfa, double *v2);

extern void VvecSparseVectors (SparseVector b1, SparseVector b2, double alfa, 
														 double *v);

extern void AddSparseVectors (SparseVector b1, SparseVector b2, 
														 ptr_SparseVector b3);

extern void JoinSparseVectors (SparseVector b1, SparseVector b2, 
															ptr_SparseVector b3);

#endif

