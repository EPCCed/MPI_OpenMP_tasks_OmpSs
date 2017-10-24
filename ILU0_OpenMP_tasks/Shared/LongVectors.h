#ifndef LongVector

#define LongVector 1

#include <stdio.h>

extern void CreateLongs (long **vlng, int dim);

extern void ReallocLongs (long **vlng, int dim);

extern void InitLongs (long *vlng, int dim, long frst, long incr);
		
extern void CopyLongs (long *src, long *dst, int dim);

extern void CopyLongsInts (long *src, int *dst, int dim);

extern void CopyIntsLongs (int *src, long *dst, int dim);

extern void CopyShiftLongs (long *src, long *dst, int dim, long shft);
	
extern void GetLongFromString (char *string, long *pnum, int dimC, int shft);

extern void GetLongsFromString (char *string, long *vlng, int dimN, int dimC, int shft);

extern void TransformLengthtoHeaderL (long *vlng, int dim);

extern void TransformHeadertoLengthL (long *vlng, int dim);

extern void ComputeHeaderfromLengthL (long *len, long *head, int dim);

extern void ComputeLengthfromHeaderL (long *head, long *len, int dim);

typedef int (*Longs_SortPerm_func) (long , long , int, int);

extern int HeapSortPermLongs_IncrPos (long n1, long n2, int p1, int p2);

extern int HeapSortPermLongs_DecrPos (long n1, long n2, int p1, int p2);

extern int HeapSortPermLongs_IncrVal (long n1, long n2, int p1, int p2);

extern int HeapSortPermLongs_DecrVal (long n1, long n2, int p1, int p2);

extern void HeapSortPermLongs_Swap (long *n1, long *n2, int *p1, int *p2);

extern void HeapSortPermLongs_Adjust (long *vlng, int *prm, int inic, int final, Longs_SortPerm_func funcSort);

extern void HeapSortPermLongs (long *vlng, int *prm, int dim, Longs_SortPerm_func funcSort);

typedef int (*Longs_Sort_func) (long , long );

extern int HeapSortLongs_IncrPos (long n1, long n2);

extern int HeapSortLongs_DecrPos (long n1, long n2);

extern void HeapSortLongs_Swap (long *n1, long *n2);

extern void HeapSortLongs_Adjust (long *vlng, int inic, int final, Longs_Sort_func funcSort);

extern void HeapSortLongs (long *vlng, int dim, Longs_Sort_func funcSort);

extern int FindSortLong (long num, long *vlng, int dim, Longs_Sort_func funcSort);
	
extern long AddLongs (long *vlng, int dim) ;

extern long AddPermuteLongs (long *vlng, int *perm, int dim);

extern long MaxLongs (long *vlng, int dim);

extern void PermuteLongs (long *vlng, int *perm, int dim);

extern void CopyPermuteLongs (long *src, int *perm, long *dst, int dim);

extern void CopyInvPermuteLongs (long *src, long *dst, int *perm, int dim);

extern void PrintLongs (long *vlng, int dim);

extern void FPrintLongs (FILE * f, long *vlng, int dim);

extern void PrintFLongs (long *vlng, int dim, int f1, int f2);

extern int ReadLongs (char *file, long **vlng);

extern void WriteLongs (char *filename, long *vlng, int dim);

extern void RemoveLongs (long **vlng);

/***********************************************************/

typedef long ** matLongs;

extern void CreateVectorPtrLongs (matLongs *mlng, int numR);

extern void CopyVectorPtrLongs (matLongs msrc, matLongs mdst, int dim);

extern void DesplVectorPtrLongs (matLongs mlng, int dim, int dspl);

extern void RemoveVectorPtrLongs (matLongs *mlng);

extern void CreateMatrixLongs (matLongs *mlng, int numR, int numC);

extern void RemoveMatrixLongs (matLongs *mlng);

#endif
