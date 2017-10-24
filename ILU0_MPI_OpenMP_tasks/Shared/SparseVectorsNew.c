#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ScalarVectors.h"
#include "SparseVectorsNew.h"

/*******************************************************************/

// This routine creates a sparse vector with nz elements at most
void CreateSparseVector (ptr_SparseVector b, int nz) {
	// The size is equal to 0, because no information is stored yet
	b->dim = 0;
	CreateInts (&(b->vpos), nz);
	CreateDoubles (&(b->vval), nz);
}

// This routine copies the data on src to dst
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
void CopySparseVectors (SparseVector src, int indexS, ptr_SparseVector dst, 
												int indexD) { 
	// Only the dimension of sparse vector is copied
	dst->dim = src.dim;
//	CopyInts (src.vpos, dst->vpos, src.dim);
	CopyShiftInts (src.vpos, dst->vpos, src.dim, indexD-indexS);
	CopyDoubles (src.vval, dst->vval, src.dim);
}

// This routine prints on the screen the data 
void PrintSparseVectors (SparseVector b) {
	int i, *pi = b.vpos;
	double *pd = b.vval;

	for (i=0; i<b.dim; i++)
		printf ("(%d,%f) ", *(pi++), *(pd++));
	printf ("\n");
}

// This routine writes on file f the data 
void FPrintSparseVectors (FILE *f, SparseVector b) {
	int i, *pi = b.vpos;
	double *pd = b.vval;

	for (i=0; i<b.dim; i++)
		fprintf (f, "(%d,%f) ", *(pi++), *(pd++));
	fprintf (f, "\n");
}

// This routine liberates the memory related to the sparse vector b
void RemoveSparseVector (ptr_SparseVector b) { 
	b->dim = -1; RemoveDoubles (&(b->vval)); RemoveInts (&(b->vpos)); 
}

/*******************************************************************/

// This function permits to define an incremental order of the index in which 
// values -1 are placed at the end
int HeapSortSparseVector_IncrPosWithNeg (int n1, int n2, double x1, double x2) { 
	return ((n2 == -1) || ((n1 != -1) && (n1 < n2))); 
}

// This function permits to define an decremental order of the index in which 
// values -1 are placed at the beginning
int HeapSortSparseVector_DecrPosWithNeg (int n1, int n2, double x1, double x2) { 
	return ((n1 == -1) || ((n2 != -1) && (n1 > n2))); 
}

// This function permits to define an incremental order of the index	
int HeapSortSparseVector_IncrPos (int n1, int n2, double x1, double x2) { 
	return (n1 < n2); 
}

// This function permits to define an decremental order of the index	
int HeapSortSparseVector_DecrPos (int n1, int n2, double x1, double x2) { 
	return (n1 > n2); 
}

// This function permits to define an incremental order of the nonzero
int HeapSortSparseVector_IncrVal (int n1, int n2, double x1, double x2) { 
	return (x1 < x2); 
}

// This function permits to define an decremental order of the nonzero
int HeapSortSparseVector_DecrVal (int n1, int n2, double x1, double x2) { 
	return (x1 > x2); 
}

// This function permits to define an incremental order of the absolute value
// of the nonzero
int HeapSortSparseVector_IncrAbsVal (int n1, int n2, double x1, double x2) { 
	return (fabs(x1) < fabs(x2)); 
}

// This function permits to define an decremental order of the absolute value
// of the nonzero
int HeapSortSparseVector_DecrAbsVal (int n1, int n2, double x1, double x2) { 
	return (fabs(x1) > fabs(x2)); 
}

// This routine permits to interchange the index and the nonzero of two 
// components of a sparse vector
void HeapSortSparseVector_Swap (int *n1, int *n2, double *x1, double *x2) {
	int n; double x;

	n = *n1; *n1 = *n2; *n2 = n;
	x = *x1; *x1 = *x2; *x2 = x;
}

// This routine is required to implement the algorithm HeapSort
void HeapSortSparseVector_Adjust (int *vpos, double *vval, int inic, int final, 
																		SparseVector_Sort_func funcSort) {
	int j, k, _bool;
		
	k = inic; j = (k << 1) + 1; _bool = 1;
	while ((j <= final) && _bool) {
		if (j < final) j += funcSort (vpos[j], vpos[j+1], vval[j], vval[j+1]);
		_bool = funcSort (vpos[k], vpos[j], vval[k], vval[j]);
		if (_bool) { 
			HeapSortSparseVector_Swap (vpos+k, vpos+j, vval+k, vval+j); 
			k = j; j = (k << 1) + 1;
		}
	}
}

// This routine implements the algorihtm HeapSort to order a sparse vector.
// The criteria is defined by the function funcSort
void HeapSortSparseVector (SparseVector b, SparseVector_Sort_func funcSort) {
	int i, ind, medio = b.dim / 2;

	for (i=1; i<=medio; i++)
		HeapSortSparseVector_Adjust (b.vpos, b.vval, (medio-i), (b.dim-1), funcSort);
	for (i=0; i<(b.dim-1); i++) {
		ind = ((b.dim-2) - i);
		HeapSortSparseVector_Swap (b.vpos, b.vpos+(ind+1), b.vval, b.vval+(ind+1));
		HeapSortSparseVector_Adjust(b.vpos, b.vval, 0, ind, funcSort);
	}
}

/*******************************************************************/

int QuickSortSparseVector_LookPivot (SparseVector b, int frs, int lst) { 
	int izq = frs+1, der = lst, *vpos = b.vpos, piv = vpos[frs]; 
	double *vval = b.vval;

	do {
		while ((vpos[izq]<=piv)&&(izq<=der)) izq++; 
		while ((vpos[der]>piv)&&(izq<=der)) der--;
		if (izq < der) {
			HeapSortSparseVector_Swap (vpos+izq, vpos+der, vval+izq, vval+der);
			izq++; der--; 
		} 
	} while (izq<=der); 
	HeapSortSparseVector_Swap (vpos+frs, vpos+der, vval+frs, vval+der);

	return der;
}

// This routine implements the algorihtm QuickSort to order a sparse vector.
// The criteria is defined by the function funcSort
void QuickSortSparseVector (SparseVector b, int frs, int lst, SparseVector_Sort_func funcSort) {
	int pos;

	if (frs < lst) {
		pos = QuickSortSparseVector_LookPivot (b, frs, lst); 
		QuickSortSparseVector (b, frs, pos-1, funcSort);
		QuickSortSparseVector (b, pos+1, lst, funcSort);
	}
}

/*******************************************************************/

// From a dense vector, this routine obtains the correponding sparse vector.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void GatherSparseVector (double *v, int n, ptr_SparseVector b, int index) {
	int i, dim = 0, *pi1 = b->vpos;
	double *pv1 = b->vval, *pv2 = v;

	for (i=0; i<n; i++) 
		if (*pv2 != 0.0) { 
			*(pi1++) = i+index; *(pv1++) = *(pv2++); dim++; 
		} else pv2++;
	b->dim = dim;
}

// From a dense vector, this routine obtains the correponding sparse vector.
// At the end, the vector is filled of zeros.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void GatherZeroSparseVector (double *v, int n, ptr_SparseVector b, int index) {
	int i, dim = 0, *pi1 = b->vpos;
	double *pv1 = b->vval, *pv2 = v;

	for (i=0; i<n; i++) 
		if (*pv2 != 0.0) { 
			*(pi1++) = i+index; *(pv1++) = *pv2; *(pv2++) = 0.0; dim++; 
		} else pv2++;
	b->dim = dim;
}

// This routine puts the information included in the sparse vector b
// on the vector v.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void ScatterSparseVector (SparseVector b, int index, double *v) {
	int i, *pi1 = b.vpos;
	double *pv1 = b.vval, *pv2 = v - index;

	for (i=0; i<b.dim; i++) pv2[*(pi1++)] = *(pv1++);
}

/*******************************************************************/

// This routine computes the dot product between the sparse vector b
// and the vector v.
// The parameter index indicates if 0-indexing or 1-indexing is used.
double DotSparseVector (SparseVector b, int index, double *v) {
	int i, *pi1 = b.vpos;
	double res = 0.0, *pv1 = b.vval, *pv2 = v - index;

	for (i=0; i<b.dim; i++) res += pv2[*(pi1++)] * *(pv1++);

	return res;
}

// This routine computes the dot product between the sparse vector b
// and the vector v. Both sparse vectors have be sorted respect the position.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
double DotSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2) {
	int i = 0, j = 0, *pi1 = b1.vpos, *pi2 = b2.vpos;
	double res = 0.0, *pv1 = b1.vval, *pv2 = b2.vval;

	while ((i < b1.dim) && (j < b2.dim))
		if ((*pi1 - index1) < (*pi2 - index2)) { 
			pi1++; pv1++; i++; 
		} else if ((*pi1 - index1) > (*pi2 - index2)) { 
			pi2++; pv2++; j++; 
		} else {
			res += *(pv1++) * *(pv2++);
			pi1++; i++; pi2++; j++;
		}

	return res;
}

// This routine computes the axpy between the sparse vector b, the scalar alpha
// and the vector v. { v += alpha * b }
// The parameter index indicates if 0-indexing or 1-indexing is used.
void AxpySparseVector (SparseVector b, int index, double alfa, double *v) {
	int i, *pi1 = b.vpos;
	double *pv1 = b.vval, *pv2 = v - index;

	for (i=0; i<b.dim; i++) pv2[*(pi1++)] += alfa * *(pv1++);
}

// This routine scales the sparse vector b by the scalar alpha
void ScalSparseVector (SparseVector b, double alfa) {
	int i = 0;
	double *pd = b.vval;

	for (i=0; i<b.dim; i++) *(pd++) *= alfa;
}

// This routine scales the sparse vector b by the scalar alpha and the vector v.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void ScalVSparseVector (double alfa, double *v, SparseVector b, int index) {
	int i = 0, *pi = b.vpos;
	double *pd = b.vval, *pv = v - index;

	for (i=0; i < b.dim; i++) { 
		*pd = alfa * *pd * pv[*pi]; pi++; pd++;
	}
}

// This routine computes the more complex axpy between the sparse vector b, the scalar alpha
// and the vectors v1 and v2. { v2 += alpha * b * v1 }
// The parameter index indicates if 0-indexing or 1-indexing is used.
void VvecSparseVector (SparseVector b, int index, double *v1, double alfa, double *v2) {
	int i = 0, *pi = b.vpos;
	double *pd = b.vval, *pv1 = v1 - index, *pv2 = v2 - index;

	for (i=0; i < b.dim; i++) { 
		pv2[*pi] += alfa * *(pd++) * pv1[*pi]; pi++; 
	}
}

// This routine computes the more complex axpy between the sparse vectors b1 and b2, 
// the scalar alpha and the vector v. { v += alpha * b1 * b2 }
// Both sparse vectors have to be sorted respect the position.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
void VvecSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2, 
												double alfa, double *v) {
	int i = 0, j = 0, *pi1 = b1.vpos, *pi2 = b2.vpos;
	double *pd1 = b1.vval, *pd2 = b2.vval, *pv = v - index1;

	while ((i < b1.dim) && (j < b2.dim))
		if ((*pi1 - index1) < (*pi2 - index2)) { 
			i++; pi1++; pd1++; 
		} else if ((*pi1 - index1) > (*pi2 - index2)) { 
			j++; pi2++; pd2++; 
		} else { 
			pv[*(pi1++)] += alfa * *(pd1++) * *(pd2++); i++; j++; pi2++; 
		}
}

// This routine computes the addition of the sparse vectors b1 and b2 on b3.
// The sparse vectors b1 and b2 have to be sorted respect the position, 
// and the sparse vector b3 is maintained ordered.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1, b2 and b3.
void AddSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2, 
												ptr_SparseVector b3, int index3) {
	int i1 = 0, i2 = 0, i3 = 0;
	int n1 = b1.dim, n2 = b2.dim;
	int index13 = index3-index1, index23 = index3-index2;
	int *pi1 = b1.vpos, *pi2 = b2.vpos, *pi3 = b3->vpos;
	double *pv1 = b1.vval, *pv2 = b2.vval, *pv3 = b3->vval;

	if ((n1 == 0) || ((pi1[n1-1] - index1) < (*pi2 - index2))) {
		// If the biggest component of b1 is smaller than the smallest component
		// of b2, both vectors are concatenated.
		CopySparseVectors (b1, index1, b3, index3); 
		b3->vpos += n1; b3->vval += n1; 
		CopySparseVectors (b2, index2, b3, index3); 
		b3->vpos -= n1; b3->vval -= n1; 
		b3->dim = (n1+n2);
	} else if ((n2 == 0) || ((pi2[n2-1] - index2) < (*pi1 - index1))) {
		// If the biggest component of b2 is smaller than the smallest component
		// of b1, both vectors are concatenated.
		CopySparseVectors (b2, index2, b3, index3); 
		b3->vpos += n2; ; b3->vval += n2;
		CopySparseVectors (b1, index1, b3, index3); 
		b3->vpos -= n2; ; b3->vval -= n2; 
		b3->dim = (n1+n2);
	} else {
		// Otherwise, both vectors have to be analyzed
		while ((i1 < n1) && (i2 < n2))
			if ((*pi1 - index1) < (*pi2 - index2)) { 
				// If the element only appears en b1, this is included.
				*(pv3++) = *(pv1++); *(pi3++) = (*(pi1++) + index13); i1++; i3++; 
			} else if ((*pi1 - index1) > (*pi2 - index2)) { 
				// If the element only appears en b2, this is included.
				*(pv3++) = *(pv2++); *(pi3++) = (*(pi2++) + index23); i2++; i3++; 
			} else {
				// If the element appears in both vectors, the addition has to be made.
				*pv3 = *(pv1++) + *(pv2++);
				if (*pv3 != 0.0) { 
					*(pi3++) = (*pi1 + index13); i3++; pv3++; 
				}
				i1++; pi1++; i2++; pi2++;
			}
		if (i1 < n1) {
			// If some elements still appear in b1, these are copied to the end
			CopyShiftInts (pi1, pi3, (n1-i1), index13);
			CopyDoubles (pv1, pv3, (n1-i1));
			i3 += (n1-i1); pi3 += (n1-i1); pv3 += (n1-i1);
		} else if (i2 < n2) {
			// If some elements still appear in b2, these are copied to the end
			CopyShiftInts (pi2, pi3, (n2-i2), index23);
			CopyDoubles (pv2, pv3, (n2-i2));
			i3 += (n2-i2); pi3 += (n2-i2); pv3 += (n2-i2);
		}
		b3->dim = i3;
	}
}

// This routine joins the sparse vectors b1 and b2 on b3, but no addition is made.
// The sparse vectors b1 and b2 have to be sorted respect the position, 
// and the sparse vector b3 is maintained ordered.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1, b2 and b3.
void JoinSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2, 
													ptr_SparseVector b3, int index3) {
	int i1 = 0, i2 = 0, i3 = 0;
	int index13 = index3-index1, index23 = index3-index2;
	int n1 = b1.dim, n2 = b2.dim;
	int *pi1 = b1.vpos, *pi2 = b2.vpos, *pi3 = b3->vpos;
	double *pv1 = b1.vval, *pv2 = b2.vval, *pv3 = b3->vval;
	
	if ((n1 == 0) || ((pi1[n1-1] - index1) < (*pi2 - index2))) {
		// If the biggest component of b1 is smaller than the smallest component
		// of b2, both vectors are concatenated.
		CopySparseVectors (b1, index1, b3, index3); 
		b3->vpos += n1; b3->vval += n1; 
		CopySparseVectors (b2, index2, b3, index3); 
		b3->vpos -= n1; b3->vval -= n1; 
		b3->dim = (n1+n2);
	} else if ((n2 == 0) || ((pi2[n2-1] - index2) < (*pi1 - index1))) {
		// If the biggest component of b2 is smaller than the smallest component
		// of b1, both vectors are concatenated.
		CopySparseVectors (b2, index2, b3, index3); 
		b3->vpos += n2; ; b3->vval += n2;
		CopySparseVectors (b1, index1, b3, index3); 
		b3->vpos -= n2; ; b3->vval -= n2; 
		b3->dim = (n1+n2);
	} else {
		// Otherwise, both vectors have to be analyzed
		while ((i1 < n1) && (i2 < n2))
			if ((*pi1 - index1) < (*pi2 - index2)) { 
				// If the element only appears en b1, this is included.
				*(pv3++) = *(pv1++); *(pi3++) = *(pi1++); i1++; i3++; 
			} else if ((*pi1 - index1) > (*pi2 - index2)) { 
				// If the element only appears en b2, this is included.
				*(pv3++) = *(pv2++); *(pi3++) = *(pi2++); i2++; i3++; 
			} else {
				// If the some position appears in both vector, 
				// the addition has to be made.
				printf ("ERRROOOOR (JoinSparseVector) \n"); exit (-1); 
			}
		if (i1 < n1) {
			// If some elements still appear in b1, these are copied to the end
			CopyShiftInts (pi1, pi3, (n1-i1), index13);
			CopyDoubles (pv1, pv3, (n1-i1));
			i3 += (n1-i1); pi3 += (n1-i1); pv3 += (n1-i1);
		} else if (i2 < n2) {
			// If some elements still appear in b2, these are copied to the end
			CopyShiftInts (pi2, pi3, (n2-i2), index23);
			CopyDoubles (pv2, pv3, (n2-i2));
			i3 += (n2-i2); pi3 += (n2-i2); pv3 += (n2-i2);
		}
		b3->dim = i3;
	}
}

// This routine computes the number of elements which appears in both sparse vectors.
// The sparse vectors b1 and b2 have to be sorted respect the position.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing or 
// 1-indexing is used to store the sparse vectors b1 and b2.
long int IntersectSparseVectors (SparseVector b1, int index1, SparseVector b2, int index2) {
	int i1 = 0, i2 = 0;
	int n1 = b1.dim, n2 = b2.dim;
	int *pi1 = b1.vpos, *pi2 = b2.vpos;
	long int Jcont = 0;
	
	if ((n1 == 0) || ((pi1[n1-1] - index1) < (*pi2 - index2))) {
		;
	} else if ((n2 == 0) || ((pi2[n2-1] - index2) < (*pi1 - index1))) {
		;
	} else {
		// Otherwise, both vectors have to be analyzed
		while ((i1 < n1) && (i2 < n2))
			if ((*pi1 - index1) < (*pi2 - index2)) { 
				i1++;
			} else if ((*pi1 - index1) > (*pi2 - index2)) { 
				i2++;
			} else {
				i1++; i2++; Jcont++;
			}
		if (i1 < n1) {
			;
		} else if (i2 < n2) {
			;
		}
	}
	return Jcont;
}

/*******************************************************************/
