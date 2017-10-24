#ifndef FloatVector

#define FloatVector 1

#include <stdio.h>

extern void CreateFloats (float **vflt, int dim);

extern void ReallocFloats (float **vftl, int dim);

extern void InitFloats (float *vflt, int dim, float frst, float incr);

extern void CopyFloats (float *src, float *dst, int dim);

extern void ScaleFloats (float *vflt, float scal, int dim);

extern void AxpyFloats (float alfa, float *vflt1, float *vflt2, int dim);

extern void XpayFloats (float *vflt1, float alfa, float *vflt2, int dim);

extern float DotFloats (float *vflt1, float *vflt2, int dim);

extern void GetFloatFromString (char *string, float *pflt, int dimC);

extern void GetFloatsFromString (char *string, float *vflt, int dimN, int dimC);

typedef int (*Floats_SortPerm_func) (float , float , int, int);

extern int HeapSortPermFloats_IncrPos (float n1, float n2, int p1, int p2);

extern int HeapSortPermFloats_DecrPos (float n1, float n2, int p1, int p2);

extern int HeapSortPermFloats_IncrVal (float n1, float n2, int p1, int p2);

extern int HeapSortPermFloats_DecrVal (float n1, float n2, int p1, int p2);

extern void HeapSortPermFloats_Swap (float *n1, float *n2, int *p1, int *p2);

extern void HeapSortPermFloats_Adjust (float *vflt, int *prm, int inic, int final, Floats_SortPerm_func funcSort);

extern void HeapSortPermFloats (float *vflt, int *prm, int dim, Floats_SortPerm_func funcSort);

typedef int (*Floats_Sort_func) (float , float );

extern int HeapSortFloats_IncrPos (float n1, float n2);

extern int HeapSortFloats_DecrPos (float n1, float n2);

extern void HeapSortFloats_Swap (float *n1, float *n2);

extern void HeapSortFloats_Adjust (float *vflt, int inic, int final, Floats_Sort_func funcSort);

extern void HeapSortFloats (float *vflt, int dim, Floats_Sort_func funcSort);

extern int FindSortFloat (float num, float *vflt, int dim, Floats_Sort_func funcSort);

extern void VdivFloats (float alfa, float *src, float *dst, int dim);

extern void VvecFloats (float alfa, float *src1, float *src2, float beta, float *dst, int dim);

extern float AddFloats (float *vflt, int dim);

extern float AddPermuteFloats (float *vflt, int *perm, int dim);

extern float MaxFloats (float *vflt, int dim);

extern void CopyPermuteFloats (float *src, int *perm, float *dst, int dim);

extern void CopyInvPermuteFloats (float *src, float *dst, int *perm, int dim);

extern void PrintFloats (float *vflt, int dim);

extern void FPrintFloats (FILE *f, float *vflt, int dim);

extern void PrintFFloats (float *vflt, int dim, int f1, int f2);

extern int ReadFloats (char *file, float **vflt);

extern void WriteFloats (char *filename, float *vflt, int dim);

extern void RemoveFloats (float **vflt);

/****************************************************************************/

typedef float ** matFloats;

extern void CreateVectorPtrFloats (matFloats *mflt, int numR);

extern void CopyVectorPtrFloats (matFloats msrc, matFloats mdst, int dim);

extern void DesplVectorPtrFloats (matFloats mflt, int dim, int dspl);

extern void RemoveVectorPtrFloats (matFloats *mflt);

extern void CreateMatrixFloats (matFloats *mflt, int numR, int numC);

extern void RemoveMatrixFloats (matFloats *mflt);

#endif

