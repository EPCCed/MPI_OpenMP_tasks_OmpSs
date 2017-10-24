#ifndef SparseVectorNewTip

#define SparseVectorNewTip 1

typedef struct SparseVector {
	int dim;
	int *vpos;
	double *vval;
	} SparseVector, *ptr_SparseVector;

/*******************************************************************/

// This routine creates a sparse vector with nz elements at most
extern void CreateSparseVector (ptr_SparseVector b, int nz);

// This routine copies the data on src to dst
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
extern void CopySparseVectors (SparseVector src, int indexS, ptr_SparseVector dst, 
												int indexD);

// This routine prints on the screen the data 
extern void PrintSparseVectors (SparseVector b);

// This routine writes on file f the data 
extern void FPrintSparseVectors (FILE *f, SparseVector b);

// This routine liberates the memory related to the sparse vector b
extern void RemoveSparseVector (ptr_SparseVector b); 

/*******************************************************************/

typedef int (*SparseVector_Sort_func) (int , int , double, double);

// This function permits to define an incremental order of the index in which 
// values -1 are placed at the end
extern int HeapSortSparseVector_IncrPosWithNeg (int n1, int n2, double x1, double x2); 

// This function permits to define an decremental order of the index in which 
// values -1 are placed at the beginning
extern int HeapSortSparseVector_DecrPosWithNeg (int n1, int n2, double x1, double x2); 

// This function permits to define an incremental order of the index	
extern int HeapSortSparseVector_IncrPos (int n1, int n2, double x1, double x2); 

// This function permits to define an decremental order of the index	
extern int HeapSortSparseVector_DecrPos (int n1, int n2, double x1, double x2); 

// This function permits to define an incremental order of the nonzero
extern int HeapSortSparseVector_IncrVal (int n1, int n2, double x1, double x2); 

// This function permits to define an decremental order of the nonzero
extern int HeapSortSparseVector_DecrVal (int n1, int n2, double x1, double x2); 

// This function permits to define an incremental order of the absolute value
// of the nonzero
extern int HeapSortSparseVector_IncrAbsVal (int n1, int n2, double x1, double x2); 

// This function permits to define an decremental order of the absolute value
// of the nonzero
extern int HeapSortSparseVector_DecrAbsVal (int n1, int n2, double x1, double x2); 

// This routine permits to interchange the index and the nonzero of two 
// components of a sparse vector
extern void HeapSortSparseVector_Swap (int *n1, int *n2, double *x1, double *x2);

// This routine is required to implement the algorithm HeapSort
extern void HeapSortSparseVector_Adjust (int *vpos, double *vval, int inic, int final, 
																					SparseVector_Sort_func funcSort);

// This routine implements the algorihtm HeapSort to order a sparse vector.
// The criteria is defined by the function funcSort
extern void HeapSortSparseVector (SparseVector b, SparseVector_Sort_func funcSort);

/*******************************************************************/

extern int QuickSortSparseVector_LookPivot (SparseVector b, int frs, int lst);

extern void QuickSortSparseVector (SparseVector b, int frs, int lst, SparseVector_Sort_func funcSort);

/*******************************************************************/

// From a dense vector, this routine obtains the correponding sparse vector.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void GatherSparseVector (double *v, int n, ptr_SparseVector b, int index);

// From a dense vector, this routine obtains the correponding sparse vector.
// At the end, the vector is filled of zeros.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void GatherZeroSparseVector (double *v, int n, ptr_SparseVector b, int index);

// This routine puts the information included in the sparse vector b
// on the vector v.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void ScatterSparseVector (SparseVector b, int index, double *v);

/*******************************************************************/

// This routine computes the dot product between the sparse vector b
// and the vector v.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern double DotSparseVector (SparseVector b, int index, double *v);

// This routine computes the dot product between the sparse vector b
// and the vector v. Both sparse vectors have be sorted respect the position.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
extern double DotSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2);

// This routine computes the axpy between the sparse vector b, the scalar alpha
// and the vector v. { v += alpha * b }
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void AxpySparseVector (SparseVector b, int index, double alfa, double *v);

// This routine scales the sparse vector b by the scalar alpha
extern void ScalSparseVector (SparseVector b, double alfa);

// This routine scales the sparse vector b by the scalar alpha and the vector v.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void ScalVSparseVector (double alfa, double *v, SparseVector b, int index);

// This routine computes the more complex axpy between the sparse vector b, the scalar alpha
// and the vectors v1 and v2. { v2 += alpha * b * v1 }
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void VvecSparseVector (SparseVector b, int index, double *v1, double alfa, double *v2);

// This routine computes the more complex axpy between the sparse vectors b1 and b2, 
// the scalar alpha and the vector v. { v += alpha * b1 * b2 }
// Both sparse vectors have to be sorted respect the position.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
extern void VvecSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2, 
																double alfa, double *v);

// This routine computes the addition of the sparse vectors b1 and b2 on b3.
// The sparse vectors b1 and b2 have to be sorted respect the position, 
// and the sparse vector b3 is maintained ordered.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1, b2 and b3.
extern void AddSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2, 
																ptr_SparseVector b3, int index3);

// This routine joins the sparse vectors b1 and b2 on b3, but no addition is made.
// The sparse vectors b1 and b2 have to be sorted respect the position, 
// and the sparse vector b3 is maintained ordered.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1, b2 and b3.
extern void JoinSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2, 
																ptr_SparseVector b3, int index3);

// This routine computes the number of elements which appears in both sparse vectors.
// The sparse vectors b1 and b2 have to be sorted respect the position.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
extern long int IntersectSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2);

/*******************************************************************/

#endif 

