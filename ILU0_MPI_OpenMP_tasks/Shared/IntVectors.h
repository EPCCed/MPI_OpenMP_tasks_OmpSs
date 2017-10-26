#ifndef IntVector

#define IntVector 1

#include <stdio.h>

extern void CreateInts (int **vint, int dim);

extern void ReallocInts (int **vint, int dim);

extern void InitInts (int *vint, int dim, int frst, int incr);

extern void GetIntFromString (char *string, int *pnum, int dimC, int shft);

extern void CopyInts (int *src, int *dst, int dim);

extern void CopyShiftInts (int *src, int *dst, int dim, int shft);

extern void CopyShiftRangsInts (int *src, int *dst, int dim1,int *rangs,  int *shfts, int dim2);

extern void GetIntsFromString (char *string, int *vint, int dimN, int dimC, int shft);

extern void GetFormatsFromString (char *string, int *vint, int dimN, int dimC);

extern void GetFormatsFromString2 (char *string, int *vint, int dimN, int dimC);

extern void TransformLengthtoHeader (int *vint, int dim);

extern void TransformHeadertoLength (int *vint, int dim);

extern void ComputeHeaderfromLength (int *len, int *head, int dim);

extern void ComputeLengthfromHeader (int *head, int *len, int dim);

typedef int (*Ints_SortPerm_func) (int , int , int, int);

extern int HeapSortPermInts_IncrPos (int n1, int n2, int p1, int p2);

extern int HeapSortPermInts_DecrPos (int n1, int n2, int p1, int p2);

extern int HeapSortPermInts_IncrVal (int n1, int n2, int p1, int p2);

extern int HeapSortPermInts_DecrVal (int n1, int n2, int p1, int p2);

extern void HeapSortPermInts_Swap (int *n1, int *n2, int *p1, int *p2);

extern void HeapSortPermInts_Adjust (int *vint, int *prm, int inic, int final, 
																		 Ints_SortPerm_func funcSort);

extern void HeapSortPermInts (int *vint, int *prm, int dim, Ints_SortPerm_func funcSort);

typedef int (*Ints_Sort_func) (int , int );

extern int HeapSortInts_IncrPos (int n1, int n2);

extern int HeapSortInts_DecrPos (int n1, int n2);

extern void HeapSortInts_Swap (int *n1, int *n2);

extern void HeapSortInts_Adjust (int *vint, int inic, int final, Ints_Sort_func funcSort);

extern void HeapSortInts (int *vint, int dim, Ints_Sort_func funcSort);

extern int FindSortInt (int num, int *vint, int dim, Ints_Sort_func funcSort);
 
extern int AddInts (int *vint, int dim); 

extern int AddPermuteInts (int *vint, int *perm, int index, int dim);

extern int MaxInts (int *vint, int dim);

extern void PermuteInts (int *vint, int *perm, int index, int dim);

extern void CopyPermuteInts (int *src, int *perm, int index, int *dst, int dim);

extern void CopyInvPermuteInts (int *src, int *dst, int *perm, int index, int dim);

extern void ComputeInvPermutation (int *perm, int *iperm, int index, int dim);

extern void PrintInts (int *vint, int dim);

extern void PrintSeqInts (int inic, int dim);

extern void FPrintInts (FILE *f, int *vint, int dim);

extern void FPrintSeqInts (FILE *f, int inic, int dim);

extern void PrintFInts (int *vint, int dim, int f1, int f2);

extern void PrintFSeqInts (int inic, int dim, int f1, int f2);

extern void FPrintFInts (FILE *f, int *vint, int dim, int f1, int f2);

extern int ReadInts (char *file, int **vint);

extern void WriteInts (char *filename, int *vint, int dim);

extern void WriteFInts (char *filename, int *vint, int dim, int f1, int f2);

extern void RemoveInts (int **vint);

/***********************************************************/

typedef int ** matInts;

extern void CreateVectorPtrInts (matInts *mint, int numR);

extern void CopyVectorPtrInts (matInts msrc, matInts mdst, int dim);

extern void DesplVectorPtrInts (matInts mint, int dim, int dspl);

extern void RemoveVectorPtrInts (matInts *mint);

extern void CreateMatrixInts (matInts *mint, int numR, int numC);

extern void RemoveMatrixInts (matInts *mint);

#endif

