#include <stdio.h>
#include <stdlib.h>
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectorsNew.h"
#include "SparseMatricesNew.h"

/*******************************************************************/

// This routine creates de sparse matrix dst from the symmetric matrix spr.
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices.
void DesymmetrizeSparseMatrices (SparseMatrix src, int indexS, ptr_SparseMatrix dst, int indexD) {
	int n = src.dim1, nnz = 0;
	int *sizes = NULL;
	int *pp1 = NULL, *pp2 = NULL, *pp3 = NULL, *pp4 = NULL, *pp5 = NULL;
	int i, j, dim, indexDS = indexD - indexS;
	double *pd3 = NULL, *pd4 = NULL;

	// The vector sizes is created and initiated
	CreateInts (&sizes, n); InitInts (sizes, n, 0, 0);
	// This loop counts the number of elements in each row
	pp1 = src.vptr; pp3 = src.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = sizes - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		// The size of the corresponding row is accumulated
		dim = (*pp2 - *pp1); pp4[i] += dim;
		// Now each component of the row is analyzed
		for (j=0; j<dim; j++) {
			// The nondiagonals elements define another element in the graph
			if (*pp3 != i) pp4[*pp3]++;
			pp3++;
		}
		pp1 = pp2++; 
	}
	
	// Compute the number of nonzeros of the new sparse matrix
	nnz = AddInts (sizes, n);
	// Create the new sparse matrix
	CreateSparseMatrix (dst, indexD, n, n, nnz, 0);
	// Fill the vector of pointers
	CopyInts (sizes, (dst->vptr) + 1, n);
	dst->vptr[0] = indexD; TransformLengthtoHeader (dst->vptr, n);
	// The vector sizes is initiated with the beginning of each row
	CopyInts (dst->vptr, sizes, n);
	// This loop fills the contents of vector vpos
	pp1 = src.vptr; pp3 = src.vpos + *pp1 - indexS; 
	pp2 = pp1 + 1 ; pp4 = dst->vpos - indexD; pp5 = sizes - indexS;
	pd3 = src.vval  + *pp1 - indexS; pd4 = dst->vval - indexD;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1);
		for (j=0; j<dim; j++) {
			// The elements in the i-th row
			pp4[pp5[i]  ] = *pp3+indexDS; 
			pd4[pp5[i]++] = *pd3; 
			if (*pp3 != i) {
				// The nondiagonals elements define another element in the graph
				pp4[pp5[*pp3]  ] = i+indexDS;
				pd4[pp5[*pp3]++] = *pd3;
			}
			pp3++; pd3++;
		}
		pp1 = pp2++;
	}
	// The memory related to the vector sizes is liberated
	RemoveInts (&sizes);
}

// This routine creates de sparse matrix dst from the symmetric matrix spr, 
// on which only the upper triangle is stored
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices.
void RestrictUpperSparseMatrices (SparseMatrix src, int indexS, ptr_SparseMatrix dst, int indexD) {
	int n = src.dim1, nnz = 0;
	int *sizes = NULL;
	int *pp1 = NULL, *pp2 = NULL, *pp3 = NULL, *pp4 = NULL, *pp5 = NULL;
	int i, j, dim, indexDS = indexD - indexS;
	double *pd3 = NULL, *pd4 = NULL;

	// The vector sizes is created and initiated
	CreateInts (&sizes, n); InitInts (sizes, n, 0, 0);
	// This loop counts the number of elements in each row
	pp1 = src.vptr; pp3 = src.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = sizes - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		// The size of the corresponding row is accumulated
		dim = (*pp2 - *pp1); pp4[i] += dim;
		// Now each component of the row is analyzed
		for (j=0; j<dim; j++) {
			// The nondiagonals elements may be changed to other row
			if (*pp3 < i) { pp4[i]--; pp4[*pp3]++; }
			pp3++;
		}
		pp1 = pp2++; 
	}
	
	// Compute the number of nonzeros of the new sparse matrix
	nnz = AddInts (sizes, n);
	// Create the new sparse matrix
	CreateSparseMatrix (dst, indexD, n, n, nnz, 0);
	// Fill the vector of pointers
	CopyInts (sizes, (dst->vptr) + 1, n);
	dst->vptr[0] = indexD; TransformLengthtoHeader (dst->vptr, n);
	// The vector sizes is initiated with the beginning of each row
	CopyInts (dst->vptr, sizes, n);
	// This loop fills the contents of vector vpos
	pp1 = src.vptr; pp3 = src.vpos + *pp1 - indexS; 
	pp2 = pp1 + 1 ; pp4 = dst->vpos - indexD; pp5 = sizes - indexS;
	pd3 = src.vval  + *pp1 - indexS; pd4 = dst->vval - indexD;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1);
		for (j=0; j<dim; j++) {
			// The nondiagonals elements may be changed to other row
			if (*pp3 >= i) {
				pp4[pp5[i]  ] = *pp3+indexDS; 
				pd4[pp5[i]++] = *pd3; 
			} else {
				pp4[pp5[*pp3]  ] = i+indexDS;
				pd4[pp5[*pp3]++] = *pd3;
			}
			pp3++; pd3++;
		}
		pp1 = pp2++;
	}
	// The memory related to the vector sizes is liberated
	RemoveInts (&sizes);
}
/*******************************************************************/

// This routine computes the product { res += spr * vec }, using sparse vectors.
// The matrix spr is symmetric, and only one of the replicated value is stored.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSymSparseMatrixVector (SparseMatrix spr, int index, double *vec, double *res) {
	int i, k;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index, *pd2 = vec, *pd3 = res;
	SparseVector baux;

	if (spr.vptr == spr.vpos) {
		// If the MSR format is used, first the diagonal has to be processed
		VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);
		for (i=index; i<spr.dim1+index; i++) {
			// A sparse vector is defined from the sparse row
			baux.dim = (*pp2-*pp1); baux.vpos = pi1; baux.vval = pd1;
			// The axpy between the column i and the i-th element of vec is computed
			AxpySparseVector (baux, index, *(pd2++), res);
			// The dot product between the row i and the vector vec is computed
			*(pd3++) += DotSparseVector (baux, index, vec);
			// Prepare the computation of the next sparse row
			pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++);
		}
	} else {
		for (i=index; i<spr.dim1+index; i++) {
			// A sparse vector is defined from the sparse row
			baux.dim = (*pp2-*pp1); baux.vpos = pi1; baux.vval = pd1;
			// Look for the diagonal element in the parse row
			k = FindSortInt (i, pi1, baux.dim, HeapSortInts_IncrPos); 
			// If exists, substract its double computation
			if (k != -1) *pd3 -= *(pd1+k) * *pd2;
			// The axpy between the column i and the i-th element of vec is computed
			AxpySparseVector (baux, index, *(pd2++), res);
			// The dot product between the row i and the vector vec is computed
			*(pd3++) += DotSparseVector (baux, index, vec);
			// Prepare the computation of the next sparse row
			pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++); 
		}
	}
}

// This routine computes the product { res += spr * vec }.
// The matrix spr is symmetric, and only one of the replicated value is stored.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSymSparseMatrixVector2 (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec - index, *pres = res - index, *pd2 = res;
	double *pd1 = spr.vval + *pp1 - index;

	if (spr.vptr == spr.vpos) {
		// If the MSR format is used, first the diagonal has to be processed
		VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);

		// After, the off-diagonals elements are processed
		for (i=index; i<spr.dim1+index; i++) {
			aux = 0.0;
			for (j=*pp1; j<*pp2; j++) {
				aux += *pd1 * pvec[*pi1];
				res[*(pi1++)] += *(pd1++) * pvec[i];
			}
			*(pd2++) += aux; pp1 = pp2++;
		}

	} else {
		for (i=index; i<spr.dim1+index; i++) {
			aux = 0.0;
			for (j=*pp1; j<*pp2; j++) {
				if (*pi1 != i)
					pres[*pi1] += *pd1 * pvec[i];
				aux += *(pd1++) * pvec[*(pi1++)];
			}
			*(pd2++) += aux; pp1 = pp2++;
		}
	}
}

// This routine computes the product { res = spr * vec }.
// The matrix spr is symmetric, and only one of the replicated value is stored.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSymSparseMatrixVector3 (SparseMatrix spr, int index, double *vec, double *res) {
	// First iniatilize the vector res to zero
	InitDoubles (res, spr.dim1, 0.0, 0.0);
	// Compute the result
	ProdSymSparseMatrixVector2 (spr, index, vec, res);
}

/*******************************************************************/

