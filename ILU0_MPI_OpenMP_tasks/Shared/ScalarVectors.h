#ifndef ScalarVector

#define ScalarVector 1

#include <stdio.h>

#include "IntVectors.h"
#include "LongVectors.h"
#include "FloatVectors.h"
#include "DoubleVectors.h"

#if 0
extern int CreateInts (int **vint, int num);

extern int InitInts (int *vint, int n, int frst, int incr);

extern void GetIntFromString (char *string, int *pnum, int numC, int shft);

extern int CopyInts (int *src, int *dst, int n);

extern int CopyShiftInts (int *src, int *dst, int n, int shft);

extern void GetIntsFromString (char *string, int *vec, int numN, int numC, int shft);

extern void GetFormatsFromString (char *string, int *vec, int numN, int numC);

extern void GetFormatsFromString2 (char *string, int *vec, int numN, int numC);

extern int TransformLengthtoHeader (int *vec, int n);

extern int TransformHeadertoLength (int *vec, int n);

extern int ComputeHeaderfromLength (int *len, int *head, int n);

extern int ComputeLengthfromHeader (int *head, int *len, int n);

typedef int (*Ints_SortPerm_func) (int , int , int, int);

extern int HeapSortPermInts_IncrPos (int n1, int n2, int p1, int p2);

extern int HeapSortPermInts_DecrPos (int n1, int n2, int p1, int p2);

extern int HeapSortPermInts_IncrVal (int n1, int n2, int p1, int p2);

extern int HeapSortPermInts_DecrVal (int n1, int n2, int p1, int p2);

extern void HeapSortPermInts_Swap (int *n1, int *n2, int *p1, int *p2);

extern void HeapSortPermInts_Adjust (int *vec, int *prm, int inic, int final, 
																		 Ints_SortPerm_func funcSort);

extern void HeapSortPermInts (int *vec, int *prm, int dim, Ints_SortPerm_func funcSort);

typedef int (*Ints_Sort_func) (int , int );

extern int HeapSortInts_IncrPos (int n1, int n2);

extern int HeapSortInts_DecrPos (int n1, int n2);

extern void HeapSortInts_Swap (int *n1, int *n2);

extern void HeapSortInts_Adjust (int *vec, int inic, int final, Ints_Sort_func funcSort);

extern void HeapSortInts (int *vec, int dim, Ints_Sort_func funcSort);

extern int FindSortInt (int num, int *vec, int dim, Ints_Sort_func funcSort);
 
extern int AddInts (int *vec, int num); 

extern int AddPermuteInts (int *vec, int *perm, int num);

extern int PermuteInts (int *vec, int *perm, int num);

extern int ComputeInvPermutation (int *perm, int *iperm, int num);

extern int PrintInts (int *vint, int num);

extern int FPrintInts (FILE *f, int *vint, int num);

extern int RemoveInts (int **vint);

extern int CreateDoubles (double **vdouble, int num);

extern int InitDoubles (double *vdouble, int n, double frst, double incr);

extern int CopyDoubles (double *src, double *dst, int n);

extern void GetDoubleFromString (char *string, double *pdbl, int numC);

extern void GetDoublesFromString (char *string, double *vec, int numN, int numC);

extern void VvecDoubles (double alfa, double *src1, double *src2, double beta, double *dst, int num);

extern int PrintDoubles (double *vdouble, int num);

extern int FPrintDoubles (FILE *f, double *vdouble, int num);

extern int RemoveDoubles (double **vdouble);

#endif

#endif

