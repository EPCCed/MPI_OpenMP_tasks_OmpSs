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
#include <stdio.h>
#include <stdlib.h>
#include "reloj.h"
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectorsNew.h"
#include "SparseMatricesNew.h"

/*******************************************************************/

// This routine creates a sparseMatrix from the next parameters
// * numR defines the number of rows
// * numC defines the number of columns
// * numE defines the number of nonzero elements
// * msr indicates if the MSR is the format used to the sparse matrix
// If msr is actived, numE doesn't include the diagonal elements
// The parameter index indicates if 0-indexing or 1-indexing is used.
void CreateSparseMatrix (ptr_SparseMatrix p_spr, int index, int numR, int numC, int numE, int msr) {
	// The scalar components of the structure are initiated
	p_spr->dim1 = numR; p_spr->dim2 = numC; 
	// Only one malloc is made for the vectors of indices
	CreateInts (&(p_spr->vptr), numE+numR+1);
	// The first component of the vectors depends on the used format
	*(p_spr->vptr) = ((msr)? (numR+1): 0) + index;
	p_spr->vpos = p_spr->vptr + ((msr)? 0: (numR+1));
	// The number of nonzero elements depends on the format used
	CreateDoubles (&(p_spr->vval), numE+(numR+1)*msr);
}

// This routine adjusts the memory size of a sparseMatrix 
void ReallocSparseMatrix (ptr_SparseMatrix p_spr) {
	int msr = (p_spr->vptr == p_spr->vpos);
	int numR = p_spr->dim1, numE = p_spr->vptr[numR]-*p_spr->vptr;
	
	// Adjust the size of the vectors of indices
	ReallocInts (&(p_spr->vptr), numE+numR+1);
	// The first component of the vectors depends on the used format
	p_spr->vpos = p_spr->vptr + ((msr)? 0: (numR+1));
	// Adjust the size of the vector of values
	ReallocDoubles (&(p_spr->vval), numE+(numR+1)*msr);
}

// This routine liberates the memory related to matrix spr
void RemoveSparseMatrix (ptr_SparseMatrix spr) {
	// First the scalar are initiated
	spr->dim1 = -1; spr->dim2 = -1; 
	// The vectors are liberated
	RemoveInts (&(spr->vptr)); RemoveDoubles (&(spr->vval)); 
}

/*******************************************************************/

// This routine prints on the screen the contents of spr in a simple way.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void PrintSparseMatrix (SparseMatrix spr, int index) {
	int i, j;

	// If the MSR format is used, first the diagonals are printed
	if (spr.vptr == spr.vpos) {
		printf ("Diagonals : \n ");
		for (i=0; i<spr.dim1; i++) printf ("%f ", spr.vval[i]); printf ("\n");
	}
	// Then, the pointers to the sparse vectors are printed
	printf ("Pointers: \n ");
	if (spr.dim1 > 0)
		for (i=0; i<=spr.dim1; i++) printf ("%d ", spr.vptr[i]); printf ("\n");
	// By rows, the pairs (index,value) are printed
	printf ("Values: \n");
	for (i=0; i<spr.dim1; i++) { 
		printf (" Row %d --> ", i+index);
		for (j=(spr.vptr[i]-index); j<(spr.vptr[i+1]-index); j++)
			printf ("(%d,%f) ", spr.vpos[j], spr.vval[j]); 
		printf ("\n");	
	}
	printf ("\n");
}

// This routine writes to the file f the contents of spr in a simple way.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void FPrintSparseMatrix (FILE *f, SparseMatrix spr, int index) {
	int i, j;

	// If the MSR format is used, first the diagonals are written
	if (spr.vptr == spr.vpos) {
		fprintf (f, "Diagonals : \n ");
		for (i=0; i<spr.dim1; i++) fprintf (f, "%f ", spr.vval[i]); fprintf (f, "\n");
	}
	// Then, the pointers to the sparse vectors are written
	fprintf (f, "Pointers: \n ");
	if (spr.dim1 > 0)
		for (i=0; i<=spr.dim1; i++) fprintf (f, "%d ", spr.vptr[i]); fprintf (f, "\n");
	// By rows, the pairs (index,value) are written
	fprintf (f, "Values: \n");
	for (i=0; i<spr.dim1; i++) { 
		fprintf (f, " Row %d --> ", i+index);
		for (j=spr.vptr[i]-index; j<spr.vptr[i+1]-index; j++)
			fprintf (f, "(%d,%f) ", spr.vpos[j], spr.vval[j]); 
		fprintf (f, "\n");	
	}
	fprintf (f, "\n");
}

// This routine prints on the screen the contents of spr, such that the data
// could be read by MATLAB.	The parameter index indicates if 0-indexing 
// or 1-indexing is used. The values fI1 and fI2 refer to the way in which
// the indices are printed, while the fD1 and fD2 are related to the values.
void PrintMatlabFSparseMatrix (SparseMatrix spr, int index, int fI1, int fI2, int fD1, int fD2) {
	char formI[10], formD[10];
	int numI = 80 / fI1, numD = 80 / fD1;
	int i, j, k;

	if (spr.dim1 > 0) {
		if (fI2 > 0) sprintf (formI, "%%%d.%dd ", fI1, fI2);
		else 				 sprintf (formI, "%%%dd ", fI1);
		if (fD2 > 0) sprintf (formD, "%%%d.%de ", fD1, fD2);
		else 				 sprintf (formD, "%%%de ", fD1);
	
		// If the MSR format is used, first the diagonals are printed
		if (spr.vptr == spr.vpos) {
			printf ("d = [ ...	\n");
			k = numD;
			for (i=0; i<spr.dim1; i++) {
				printf (formD, spr.vval[i]); 
				if ((--k) == 0) { printf (" ... \n"); k = numD; }
			}
			printf ("];\n");
		}
		// Then, the pointers to the sparse vectors are printed
		printf ("vptr = [ ...	\n");
		k = numI;
		for (i=0; i<=spr.dim1; i++) {
			printf (formI, spr.vptr[i]); 
			if ((--k) == 0) { printf (" ... \n"); k = numI; }
		}
		printf ("] + %d;\n", (1-index));
		// After, the rows of the nonzeos are printed
		printf ("rows = [ ...	\n");
		k = numI;
		for (i=0; i<spr.dim1; i++) {
			for (j = spr.vptr[i]; j<spr.vptr[i+1]; j++) {
				printf (formI, i+index); 
				if ((--k) == 0) { printf (" ... \n"); k = numI; }
			}
		}
		printf ("] + %d;\n", (1-index));
		// Before, the columns of the nonzeros are printed
		printf ("cols = [ ...	\n");
		k = numI;
		for (i=0; i<spr.vptr[spr.dim1]-index; i++) {
			printf (formI, spr.vpos[i]+index); 
			if ((--k) == 0) { printf (" ... \n"); k = numI; }
		}
		printf ("] + %d;\n", (1-index));
		// Finally, the values of the nonzeros are printed
		printf ("vals = [ ...	\n");
		k = numD;
		for (i=0; i<spr.vptr[spr.dim1]-index; i++) {
			printf (formD, spr.vval[i]); 
			if ((--k) == 0) { printf (" ... \n"); k = numD; }
		}
		printf ("];\n");
		printf ("SpMat = sparse (rows, cols, vals);\n");
		if (spr.vptr == spr.vpos) {
			printf ("SpMat = SpMat + diag (d);\n");
		}
	}
}

// This routine writes on the file f the contents of spr, such that the data
// could be read by MATLAB.	The parameter index indicates if 0-indexing 
// or 1-indexing is used. The values fI1 and fI2 refer to the way in which
// the indices are printed, while the fD1 and fD2 are related to the values.
void FPrintMatlabFSparseMatrix (FILE *f, SparseMatrix spr, int index, 
																int fI1, int fI2, int fD1, int fD2) {
	char formI[10], formD[10];
	int numI = 80 / fI1, numD = 80 / fD1;
	int i, j, k;

	if (spr.dim1 > 0) {
		if (fI2 > 0) sprintf (formI, "%%%d.%dd ", fI1, fI2);
		else 				 sprintf (formI, "%%%dd ", fI1);
		if (fD2 > 0) sprintf (formD, "%%%d.%de ", fD1, fD2);
		else 				 sprintf (formD, "%%%de ", fD1);
	
		// If the MSR format is used, first the diagonals are printed
		if (spr.vptr == spr.vpos) {
			fprintf (f, "d = zeros(%d,1);	\n", spr.dim1);
			fprintf (f, "d = [ ...	\n");
			k = numD;
			for (i=0; i<spr.dim1; i++) {
				fprintf (f, formD, spr.vval[i]); 
				if ((--k) == 0) { fprintf (f, " ... \n"); k = numD; }
			}
			fprintf (f, "];\n");
		}
		// Then, the pointers to the sparse vectors are printed
		fprintf (f, "vptr = zeros(%d,1);	\n", spr.dim1+1);
		fprintf (f, "vptr = [ ...	\n");
		k = numI;
		for (i=0; i<=spr.dim1; i++) {
			fprintf (f, formI, spr.vptr[i]); 
			if ((--k) == 0) { fprintf (f, " ... \n"); k = numI; }
		}
		fprintf (f, "] + %d;\n", (1-index));
		// After, the rows of the nonzeos are printed
		fprintf (f, "rows = zeros(%d,1);	\n", spr.vptr[spr.dim1]-index);
		fprintf (f, "rows = [ ...	\n");
		k = numI;
		for (i=0; i<spr.dim1; i++) {
			for (j = spr.vptr[i]; j<spr.vptr[i+1]; j++) {
				fprintf (f, formI, i+index); 
				if ((--k) == 0) { fprintf (f, " ... \n"); k = numI; }
			}
		}
		fprintf (f, "] + %d;\n", (1-index));
		// Before, the columns of the nonzeros are printed
		fprintf (f, "cols = zeros(%d,1);	\n", spr.vptr[spr.dim1]-index);
		fprintf (f, "cols = [ ...	\n");
		k = numI;
		for (i=0; i<spr.vptr[spr.dim1]-index; i++) {
			fprintf (f, formI, spr.vpos[i]); 
			if ((--k) == 0) { fprintf (f, " ... \n"); k = numI; }
		}
		fprintf (f, "] + %d;\n", (1-index));
		// Finally, the values of the nonzeros are printed
		fprintf (f, "vals = zeros(%d,1);	\n", spr.vptr[spr.dim1]-index);
		fprintf (f, "vals = [ ...	\n");
		k = numD;
		for (i=0; i<spr.vptr[spr.dim1]-index; i++) {
			fprintf (f, formD, spr.vval[i]); 
			if ((--k) == 0) { fprintf (f, " ... \n"); k = numD; }
		}
		fprintf (f, "];\n");
		fprintf (f, "SpMat = sparse (rows, cols, vals);\n");
		if (spr.vptr == spr.vpos) {
			fprintf (f, "SpMat = SpMat + diag (d);\n");
		}
	}
}

// This routine writes on the file filename the contents of spr, such that the data
// could be read by MATLAB.	The parameter index indicates if 0-indexing 
// or 1-indexing is used. The values fI1 and fI2 refer to the way in which
// the indices are printed, while the fD1 and fD2 are related to the values.
void WriteMatlabFSparseMatrix (char *filename, SparseMatrix spr, int index, 
																int fI1, int fI2, int fD1, int fD2) {
  FILE *f = NULL;

  printf ("Writing %s of size (%d,%d,%d)\n", filename, 
							spr.dim1, spr.dim2, spr.vptr[spr.dim1]-index);
	f = OpenFile (filename, "w");
	FPrintMatlabFSparseMatrix (f, spr, index, fI1, fI2, fD1, fD2);
  fclose (f);
  printf ("Written %s\n", filename);

}

// This routine writes different files to store the contents of spr.
// The name of these files begins with the string prefix.
// The parameter index indicates if 0-indexing or 1-indexing is used. 
// The values fI1 and fI2 refer to the way in which the indices are printed, 
// while the fD1 and fD2 are related to the values.
void WriteMatlabFSparseMatrix2 (char *prefix, SparseMatrix spr, int index, 
																int fI1, int fI2, int fD1, int fD2) {
	char filename[80];

  printf ("Writing %s of size (%d,%d,%d)\n", prefix, 
							spr.dim1, spr.dim2, spr.vptr[spr.dim1]-index);

	sprintf (filename, "%s_vptr.txt", prefix);
	WriteFInts (filename, spr.vptr, spr.dim1+1, fI1, fI2);
	sprintf (filename, "%s_vpos.txt", prefix);
	WriteFInts (filename, spr.vpos, spr.vptr[spr.dim1]-index, fI1, fI2);
	sprintf (filename, "%s_vval.txt", prefix);
	WriteFDoubles (filename, spr.vval, spr.vptr[spr.dim1]-index, fD1, fD2);

  printf ("Written %s\n", prefix);

}

/*******************************************************************/

// If it is necessary, this routine changes the indexing of spr,
// from 0-indexing to 1-indexing or in the reverse case.
// The parameter lngs defines the change, where its first bit
// defines the indexing at the input and its second bit defines
// the deseared indexing at the output, thus
// * if is equal to 1, changes to 1 to 0 indexing.
// * if is equal to 2, changes to 0 to 1 indexing.
void ConvertSparseMatrix (SparseMatrix spr, int lngs) {
	int i;

	if ((spr.dim1 > 0) && (spr.dim2 > 0))
		// If the matrix is not empty
		if (((lngs & 2)>>1) != (lngs & 1)) {
			// If the indexing has to be changed
			if (lngs & 1) {
				// The change is from 1-indexing to 0-indexing
				printf ("From 1-indexing to 0-indexing\n");
				for (i=0; i<=spr.dim1; i++) spr.vptr[i]--;
				for (i=spr.vptr[0]; i<spr.vptr[spr.dim1]; i++) spr.vpos[i]--;
			} else {
				// The change is from 0-indexing to 1-indexing
				printf ("From 0-indexing to 1-indexing\n");
				for (i=spr.vptr[0]; i<spr.vptr[spr.dim1]; i++) spr.vpos[i]++;
				for (i=0; i<=spr.dim1; i++) spr.vptr[i]++;
			}
		}
}

/*******************************************************************/

// This routine sorts the rows of spr, according to the criteria defined in funcSort.
// This operation is made row by row.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void HeapSortSparseMatrix (SparseMatrix spr, int index, SparseVector_Sort_func funcSort) {
	int i; 
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;
	SparseVector baux;

	for (i=0; i<spr.dim1; i++) {
		// Defines the sparse vector to be sorted
		baux.dim = *pp2-*pp1; baux.vpos = pi1; baux.vval = pd1;
		HeapSortSparseVector (baux, funcSort);
		// Determines the values of the next iteration
		pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++); 
	}
}

// This routine permuts the columns of spr, if perm is not null,
// or it only sort the rows according to the criteria defined in funcSort.
// The parameter index indicates if 0-indexing or 1-indexing is used,
// for both, the matrix spy and the vector perm.
void PermuteColsSparseMatrix (SparseMatrix spr, int index, int *perm) {
	int i; 
	int *pp1 = spr.vptr, *pp2 = pp1 + 1, *pi1 = spr.vpos + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;
	SparseVector baux;
	double te = 0.0, te1, te2, tu = 0.0, tu1, tu2;

	// Apply a permutation on the columns if the vector perm exists
	if (perm != NULL) PermuteInts (pi1, perm, index, *(pp1+spr.dim1)-*pp1); 

	reloj (&te1, &tu1);
	for (i=0; i<spr.dim1; i++) {
		// Defines the sparse vector to be sorted
		baux.dim = *pp2-*pp1; baux.vpos = pi1; baux.vval = pd1;
		HeapSortSparseVector (baux, HeapSortSparseVector_IncrPos);
		// Determines the values of the next iteration
		pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++);
	}
	reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);
	printf ("(B) dim = %d, nnz = %6d , te = %20.10f, tu = %20.10f\n", spr.dim1, spr.vptr[spr.dim1], te, tu);
}

// This routine permuts the columns of spr, if perm is not null,
// or it only sort the rows according to the criteria defined in funcSort.
// The perm contains negative indices to mark the columns to be removed.
// The parameter index indicates if 0-indexing or 1-indexing is used,
// for both, the matrix spy and the vector perm.
void PermuteColsWithNegSparseMatrix (SparseMatrix spr, int index, int *perm) {
	int i, k1 = 0, k2 = 0; 
	int *pp1 = spr.vptr, *pp2 = pp1 + 1, *pi1 = spr.vpos + *pp1 - index, *pi2 = NULL;
	double *pd1 = spr.vval + *pp1 - index;
	SparseVector baux;

	// Apply a permutation on the columns if the vector perm exists
	if (perm != NULL) PermuteInts (pi1, perm, index, *(pp1+spr.dim1)-*pp1); 

	for (i=0; i<spr.dim1; i++) {
		// Defines the sparse vector to be sorted
		baux.dim = *pp2-*pp1; baux.vpos = pi1; baux.vval = pd1;
		HeapSortSparseVector (baux, HeapSortSparseVector_IncrPosWithNeg);
		// Locate the number of negative indices, k1
		k1 = baux.dim; pi2 = pi1+baux.dim; while ((k1 > 0) && (*(--pi2) < 0)) k1--;
		if (k2 > 0) {
			// Remove the nondesired elements of the data vectors
			CopyInts (pi1, pi1-k2, k1);
			CopyDoubles (pd1, pd1-k2, k1);
			*pp1 -= k2;
		}
		// Determines the values of the next iteration
		k2 += (baux.dim - k1); pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++);
	}
	*pp1 -= k2;
}
 
/*******************************************************************/

// This routine extracts the diagonal of spr on the vector diag.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void GetDiagonalSparseMatrix (SparseMatrix spr, int index, double *diag) {
	int i, j, dim = (spr.dim1 < spr.dim2)? spr.dim1 : spr.dim2;
	int *pp1 = NULL, *pp2 = NULL, *pi1 = NULL; 
	double *pd1 = NULL, *pd2 = diag;

	if (spr.vptr == spr.vpos)
		// If MSR format is used, only a copy has to be made
		CopyDoubles (spr.vval, diag, spr.dim1);
	else {
		// Determines the auxiliar vectors
		pp1 = spr.vptr; pp2 = pp1 + 1; 
		pi1 = spr.vpos + *pp1 - index; pd1 = spr.vval + *pp1 - index; 
		for (i=index; i<(dim+index); i++) {
			// Looks for the diagonal value
			j = (*pp2-*pp1); while ((j > 0) && (*pi1 < i)) { pi1++; pd1++; j--; }
			// Fixes the diagonal, putting 0 if it doesn't exist
			*(pd2++) = ((j > 0) && (*pi1 == i))? *pd1: 0.0;
			// Determines the values of the next iteration
			pi1 += j; pd1 += j; pp1 = (pp2++); 
		}
	}
}

// This routine extracts the diagonal of spr on the vector diag.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void GetDiagonalSparseMatrixDspls (SparseMatrix spr, int index, double *diag, int dspls) {
	int i, j, dim = (spr.dim1 < spr.dim2)? spr.dim1 : spr.dim2;
	int *pp1 = NULL, *pp2 = NULL, *pi1 = NULL; 
	double *pd1 = NULL, *pd2 = diag;

	if (spr.vptr == spr.vpos)
		// If MSR format is used, only a copy has to be made
		CopyDoubles (spr.vval+dspls, diag, spr.dim1-dspls);
	else {
		// Determines the auxiliar vectors
		pp1 = spr.vptr+dspls; pp2 = pp1+1;
		pi1 = spr.vpos + *pp1 - index; pd1 = spr.vval + *pp1 - index; 
		for (i=dspls+index; i<dim+index; i++) {
			// Looks for the diagonal value
			j = (*pp2-*pp1); while ((j > 0) && (*pi1 < i)) { pi1++; pd1++; j--; }
			// Fixes the diagonal, putting 0 if it doesn't exist
			*(pd2++) = ((j > 0) && (*pi1 == i))? *pd1: 0.0;
			// Determines the values of the next iteration
			pi1 += j; pd1 += j; pp1 = (pp2++);
		}
	}
}

// This routine computes the graph related to the symmetric sparse matrix spr, in which 
// the diagonals are not considered, and the result is stored on the vectors vptr and vpos.
// The parameters indexS and indexG indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrix and the graph.
void GetGraphSparseMatrix (SparseMatrix spr, int indexS, int *vptr, int *vpos, int indexG) {
	int n = spr.dim1;
	int *ptrSpr = NULL;
	int *pp1 = NULL, *pp2 = NULL, *pp3 = NULL, *pp4 = NULL, *pp5 = NULL;
	int i, j, dim, indexGS = indexG - indexS;

	// The vptr is initiated
	*vptr = indexG; InitInts (vptr+1, n, 0, 0);
	// This loop counts the number of elements in each row
	pp1 = spr.vptr; pp3 = spr.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = vptr + 1 - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		// The size of the corresponding row is accumulated
		dim = (*pp2 - *pp1); pp4[i] += dim;
		// Now each component of the row is analyzed
		for (j=0; j<dim; j++) {
			// The diagonal isn't considered.
			if (*pp3 == i) pp4[i]--; 
			// Otherwise, the symmetric element has to be considered
			else        pp4[*pp3]++;
			pp3++;
		}
		pp1 = pp2++; 
	}
	// Knowing the lengths, it is easy to obtain the vector vptr
	TransformLengthtoHeader (vptr, n);

	// An auxiliar vector is created to traverse spr
	CreateInts (&ptrSpr, n+1);
	// The auxiliar vector is initiated with the beginning of each row
	CopyInts (vptr, ptrSpr, n+1);
	// This loop fills the contents of vector vpos
	pp1 = spr.vptr; pp3 = spr.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = vpos - indexG; pp5 = ptrSpr - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1);
		for (j=0; j<dim; j++) {
			// The diagonal isn't considered.
			if (*pp3 != i) {
				// The nondiagonal elements define two positions in the graph
				pp4[pp5[i]++] = *pp3+indexGS; 
				pp4[pp5[*pp3]++] = i+indexGS;
			}
			pp3++;
		}
		pp1 = pp2++;
	}
	// The memory related to the auxiliar vector is liberated
	RemoveInts(&ptrSpr);
}

// This routine computes the graph related to the symmetric sparse matrix spr, in which 
// the diagonals are considered, and the result is stored on the vectors vptr and vpos.
// The parameters indexS and indexG indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrix and the graph.
void GetGraphWithDiagSparseMatrix (SparseMatrix spr, int indexS, int *vptr, int *vpos, int indexG) {
	int n = spr.dim1;
	int *ptrSpr = NULL;
	int *pp1 = NULL, *pp2 = NULL, *pp3 = NULL, *pp4 = NULL, *pp5 = NULL;
	int i, j, dim, indexGS = indexG - indexS;

	// The vptr is initiated
	*vptr = indexG; InitInts (vptr+1, n, 0, 0);
	// This loop counts the number of elements in each row
	pp1 = spr.vptr; pp3 = spr.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = vptr + 1 - indexS;
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
	// Knowing the lengths, it is easy to obtain the vector vptr
	TransformLengthtoHeader (vptr, n);

	// An auxiliar vector is created to traverse spr
	CreateInts (&ptrSpr, n+1);
	// The auxiliar vector is initiated with the beginning of each row
	CopyInts (vptr, ptrSpr, n+1);
	// This loop fills the contents of vector vpos
	pp1 = spr.vptr; pp3 = spr.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = vpos - indexG; pp5 = ptrSpr - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1);
		for (j=0; j<dim; j++) {
			// The elements in the i-th row
			pp4[pp5[i]++] = *pp3+indexGS; 
			if (*pp3 != i) {
				// The nondiagonals elements define another element in the graph
				pp4[pp5[*pp3]++] = i+indexGS;
			}
			pp3++;
		}
		pp1 = pp2++;
	}
	// The memory related to the auxiliar vector is liberated
	RemoveInts(&ptrSpr);
}

/*******************************************************************/

// This routine computes the diagonal scaling { A = diag(rowscal)*A*diag(colscal) }, 
// using sparse vectors.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ScalSparseMatrix (SparseMatrix spr, int index, double *rowscal, double *colscal) {
	int i;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index, *pd2 = rowscal;
	SparseVector baux;

	// If the MSR format is used, first the diagonal has to be processed
	if (spr.vptr == spr.vpos) {
		VvecDoubles (1.0, rowscal, spr.vval, 0.0, spr.vval, spr.dim1);
		VvecDoubles (1.0, colscal, spr.vval, 0.0, spr.vval, spr.dim1);
	}

	for (i=0; i<spr.dim1; i++) {
		// A sparse vector is defined from the sparse row on which
    // the column scale is applied. Aditionally, each sparse row
    // has to be scaled by the correspondent component of row scale.
		// compute its dot product and the vector vec
		baux.dim = (*pp2-*pp1); baux.vpos = pi1; baux.vval = pd1;
		ScalVSparseVector (*(pd2++), colscal, baux, index);
		// Prepare the computation of the next sparse row
		pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++); 
	}
}

// This routine computes the product { res += spr * vec }, using sparse vectors.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVector (SparseMatrix spr, int index, double *vec, double *res) {
	int i;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index, *pd2 = res;
	SparseVector baux;

	// If the MSR format is used, first the diagonal has to be processed
	if (spr.vptr == spr.vpos)
		VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);

	for (i=0; i<spr.dim1; i++) {
		// A sparse vector is defined from the sparse row on which
    // the dot product is computed with the vector vec.
		baux.dim = (*pp2-*pp1); baux.vpos = pi1; baux.vval = pd1;
		*(pd2++) += DotSparseVector (baux, index, vec);
		// Prepare the computation of the next sparse row
		pi1 += baux.dim; pd1 += baux.dim; pp1 = (pp2++); 
	}
}

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVector2 (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec - index, *pd2 = res;
	double *pd1 = spr.vval + *pp1 - index;

	// If the MSR format is used, first the diagonal has to be processed
	if (spr.vptr == spr.vpos)
		VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);

	for (i=0; i<spr.dim1; i++) {
		// The dot product between the row i and the vector vec is computed
		aux = 0.0;
		for (j=*pp1; j<*pp2; j++)
			aux += *(pd1++) * pvec[*(pi1++)];
		// Accumulate the obtained value on the result
		*(pd2++) += aux; pp1 = pp2++;
	}
}

// This routine computes the product { res = spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVector3 (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec - index, *pd2 = res;
	double *pd1 = spr.vval + *pp1 - index;

	// If the MSR format is used, first the diagonal has to be processed
	if (spr.vptr == spr.vpos)
		VvecDoubles (1.0, spr.vval, vec, 0.0, res, spr.dim1);

	for (i=0; i<spr.dim1; i++) {
		// The dot product between the row i and the vector vec is computed
		aux = 0.0;
		for (j=*pp1; j<*pp2; j++)
			aux += *(pd1++) * pvec[*(pi1++)];
		// Assign the obtained value on the result
		*(pd2++) = aux; pp1 = pp2++;
	}
}

/*******************************************************************/

// This routine computes the addition { dst = spr1(perm1,perm1) + spr2(perm2,perm2) }
// The vector iperm1 contains the inverse permutation of perm1, in which the value
// -1 indicates that the correspondient position doesn't exist in perm1.
// The same for iperm2 and perm2.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing 
// or 1-indexing is used to store the different sparse matrices.
void AddSparseMatricesNew (SparseMatrix src1, int index1, int *perm1, int *iperm1,
												SparseMatrix src2, int index2, int *perm2, int *iperm2,
												ptr_SparseMatrix dst, int index3) {
	int i, j, k;
	int *vlen = dst->vptr+1, *vpos = dst->vpos + *(dst->vptr) - index3;
	double diag, *vval = dst->vval + *(dst->vptr) - index3;
	SparseVector b1, b2, b;
	int    *pp1 = src1.vptr - index1, *pp2 = src2.vptr - index2;
	int    *pi1 = src1.vpos - index1, *pi2 = src2.vpos - index2;
	double *pd1 = src1.vval - index1, *pd2 = src2.vval - index2;
	double te = 0.0, te1, te2, tu = 0.0, tu1, tu2;

	// Apply the permutations required to make the addition
	reloj (&te1, &tu1);
	reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

	printf ("index1 = %d , index2 = %d , index3 = %d\n", index1, index2, index3);

	for (i=0; i<dst->dim1; i++) {
		// Initialization of the elements on which the addition will be made
		diag = 0.0; b.dim = 0; b.vpos = vpos; b.vval = vval;
		if ((iperm1[i] != -1) && (iperm2[i] != -1)) {
			// If the row appears in both matrices, the addition has to be made
			// If the MSR is used, the diag is changed
			if (src1.vptr == src1.vpos) diag += pd1[iperm1[i]];
			if (src2.vptr == src2.vpos) diag += pd2[iperm2[i]];
			// The sparse vector b1 is made from the i-th row of src1
			b1.dim = pp1[iperm1[i]+1] - pp1[iperm1[i]];
			b1.vpos = pi1 + pp1[iperm1[i]]; b1.vval = pd1 + pp1[iperm1[i]];

			PermuteInts (b1.vpos, perm1, index1, b1.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b1, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

			// The sparse vector b2 is made from the i-th row of src2
			b2.dim = pp2[iperm2[i]+1] - pp2[iperm2[i]];
			b2.vpos = pi2 + pp2[iperm2[i]]; b2.vval = pd2 + pp2[iperm2[i]];

			PermuteInts (b2.vpos, perm2, index2, b2.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b2, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

			// The addition is computed
			AddSparseVectors (b1, index1, b2, index2, &b, index3);

			PermuteInts (b1.vpos, iperm1, index1, b1.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b1, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

			PermuteInts (b2.vpos, iperm2, index2, b2.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b2, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

		} else if (iperm1[i] != -1) {
			// If the row only in matrix src1, it is copied
			// If the MSR is used, the diag is changed
			if (src1.vptr == src1.vpos) diag += pd1[iperm1[i]];
			// The sparse vector b1 is made from the i-th row of src1
			b1.dim = pp1[iperm1[i]+1] - pp1[iperm1[i]];
			b1.vpos = pi1 + pp1[iperm1[i]]; b1.vval = pd1 + pp1[iperm1[i]];
			// The copy is made

			PermuteInts (b1.vpos, perm1, index1, b1.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b1, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

			CopySparseVectors (b1, index1, &b, index3);

			PermuteInts (b1.vpos, iperm1, index1, b1.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b1, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

		} else { 
			// If the row only in matrix src2, it is copied
			// If the MSR is used, the diag is changed
			if (src2.vptr == src2.vpos) diag += pd2[iperm2[i]];
			// The sparse vector b2 is made from the i-th row of src2
			b2.dim = pp2[iperm2[i]+1] - pp2[iperm2[i]];
			b2.vpos = pi2 + pp2[iperm2[i]]; b2.vval = pd2 + pp2[iperm2[i]];
			// The copy is made

			PermuteInts (b2.vpos, perm2, index2, b2.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b2, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

			CopySparseVectors (b2, index2, &b, index3);

			PermuteInts (b2.vpos, iperm2, index2, b2.dim);
			reloj (&te1, &tu1);
			HeapSortSparseVector (b2, HeapSortSparseVector_IncrPos);
			reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

		}

		if ((src1.vptr != src1.vpos) || (src2.vptr != src2.vpos)) {
			// If some of the source matrices don't use the MSR format
			// the position of the diag has to be found
			k = FindSortInt (i+index3, vpos, b.dim, HeapSortInts_IncrPos); 
			if (k != -1) {
				// if the diagonal has been found
				if (dst->vptr == dst->vpos)	{
					// If the solution uses the MSR format
					// to accumulate the value on diag
					*(dst->vval+i) = diag + *(vval+k); b.dim--;
					// And remove the element from the sparse vector
					for (j=k; j<b.dim; j++) { 
						vpos[j] = vpos[j+1]; vval[j] = vval[j+1]; 
					}
				} else if (diag != 0.0)
					// Otherwise, if the value of diag is not null, accumulate it
					*(vval+k) += diag;
			} else if (dst->vptr == dst->vpos) 
				// If the solution uses the MSR format to assign the value on diag
				*(dst->vval+i) = diag;
			else if (diag != 0.0) {
				// Otherwise, if the value of diag is not null, copy it in the
				// correct position to mantain the order
				j = b.dim;
				while ((j > 0) && (vpos[j-1] > (i+index3))) { 
					vpos[j] = vpos[j-1]; vval[j] = vval[j-1]; j--; 
				}
				vpos[j] = i+index3; vval[j] = diag; b.dim++;
			}
		} else if (dst->vptr == dst->vpos)	
			// If the solution uses the MSR format to assign the value on diag
			*(dst->vval+i) = diag;
		else if (diag != 0.0) {
			// Otherwise, if the value of diag is not null, copy it in the
			// correct position to mantain the order
			j = b.dim;
			while ((j > 0) && (vpos[j-1] > (i+index3))) { 
				vpos[j] = vpos[j-1]; vval[j] = vval[j-1]; j--; 
			}
			vpos[j] = (i+index3); vval[j] = diag; b.dim++;
		}
		// Update the values to the next iteration
		vpos += b.dim; vval += b.dim; *(vlen++) = b.dim;
	}
	// Knowing the lengths, it is easy to obtain the vector vptr
	TransformLengthtoHeader (dst->vptr, dst->dim1);
	// Remove the permutations applied to make the addition
	reloj (&te1, &tu1);
//	PermuteColsSparseMatrix (src2, index2, iperm2); PermuteColsSparseMatrix (src1, index1, iperm1);
	reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);
}

/*******************************************************************/

// This routine computes the addition { dst = spr1(perm1,perm1) + spr2(perm2,perm2) }
// The vector iperm1 contains the inverse permutation of perm1, in which the value
// -1 indicates that the correspondient position doesn't exist in perm1.
// The same for iperm2 and perm2.
// The parameters index1, index2 and index3 indicate, respectivaly, if 0-indexing 
// or 1-indexing is used to store the different sparse matrices.
// permCols determinates if the vpos vectors of the matrices have to be permuted.
void AddSparseMatrices (SparseMatrix src1, int index1, int *perm1, int *iperm1,
												SparseMatrix src2, int index2, int *perm2, int *iperm2,
												ptr_SparseMatrix dst, int index3, int perm_cols) {
	int i, j, k;
	int *vlen = dst->vptr+1, *vpos = dst->vpos + *(dst->vptr) - index3;
	double diag, *vval = dst->vval + *(dst->vptr) - index3;
	SparseVector b1, b2, b;
	int    *pp1 = src1.vptr - index1, *pp2 = src2.vptr - index2;
	int    *pi1 = src1.vpos - index1, *pi2 = src2.vpos - index2;
	double *pd1 = src1.vval - index1, *pd2 = src2.vval - index2;
	double te = 0.0, te1, te2, tu = 0.0, tu1, tu2;

	// Apply the permutations required to make the addition
	reloj (&te1, &tu1);
	if (perm_cols) {
		PermuteColsSparseMatrix (src1, index1, perm1); PermuteColsSparseMatrix (src2, index2, perm2);
	}
	reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);

	for (i=0; i<dst->dim1; i++) {
		// Initialization of the elements on which the addition will be made
		diag = 0.0; b.dim = 0; b.vpos = vpos; b.vval = vval;
		if ((iperm1[i] != -1) && (iperm2[i] != -1)) {
			// If the row appears in both matrices, the addition has to be made
			// If the MSR is used, the diag is changed
			if (src1.vptr == src1.vpos) diag += pd1[iperm1[i]];
			if (src2.vptr == src2.vpos) diag += pd2[iperm2[i]];
			// The sparse vector b1 is made from the i-th row of src1
			b1.dim = pp1[iperm1[i]+1] - pp1[iperm1[i]];
			b1.vpos = pi1 + pp1[iperm1[i]]; b1.vval = pd1 + pp1[iperm1[i]];
			// The sparse vector b2 is made from the i-th row of src2
			b2.dim = pp2[iperm2[i]+1] - pp2[iperm2[i]];
			b2.vpos = pi2 + pp2[iperm2[i]]; b2.vval = pd2 + pp2[iperm2[i]];
			// The addition is computed
			AddSparseVectors (b1, index1, b2, index2, &b, index3);
		} else if (iperm1[i] != -1) {
			// If the row only in matrix src1, it is copied
			// If the MSR is used, the diag is changed
			if (src1.vptr == src1.vpos) diag += pd1[iperm1[i]];
			// The sparse vector b1 is made from the i-th row of src1
			b1.dim = pp1[iperm1[i]+1] - pp1[iperm1[i]];
			b1.vpos = pi1 + pp1[iperm1[i]]; b1.vval = pd1 + pp1[iperm1[i]];
			// The copy is made
			CopySparseVectors (b1, index1, &b, index3);
		} else { 
			// If the row only in matrix src2, it is copied
			// If the MSR is used, the diag is changed
			if (src2.vptr == src2.vpos) diag += pd2[iperm2[i]];
			// The sparse vector b2 is made from the i-th row of src2
			b2.dim = pp2[iperm2[i]+1] - pp2[iperm2[i]];
			b2.vpos = pi2 + pp2[iperm2[i]]; b2.vval = pd2 + pp2[iperm2[i]];
			// The copy is made
			CopySparseVectors (b2, index2, &b, index3);
		}

		if ((src1.vptr != src1.vpos) || (src2.vptr != src2.vpos)) {
			// If some of the source matrices don't use the MSR format
			// the position of the diag has to be found
			k = FindSortInt (i+index3, vpos, b.dim, HeapSortInts_IncrPos); 
			if (k != -1) {
				// if the diagonal has been found
				if (dst->vptr == dst->vpos)	{
					// If the solution uses the MSR format
					// to accumulate the value on diag
					*(dst->vval+i) = diag + *(vval+k); b.dim--;
					// And remove the element from the sparse vector
					for (j=k; j<b.dim; j++) { 
						vpos[j] = vpos[j+1]; vval[j] = vval[j+1]; 
					}
				} else if (diag != 0.0)
					// Otherwise, if the value of diag is not null, accumulate it
					*(vval+k) += diag;
			} else if (dst->vptr == dst->vpos) 
				// If the solution uses the MSR format to assign the value on diag
				*(dst->vval+i) = diag;
			else if (diag != 0.0) {
				// Otherwise, if the value of diag is not null, copy it in the
				// correct position to mantain the order
				j = b.dim;
				while ((j > 0) && (vpos[j-1] > (i+index3))) { 
					vpos[j] = vpos[j-1]; vval[j] = vval[j-1]; j--; 
				}
				vpos[j] = i+index3; vval[j] = diag; b.dim++;
			}
		} else if (dst->vptr == dst->vpos)	
			// If the solution uses the MSR format to assign the value on diag
			*(dst->vval+i) = diag;
		else if (diag != 0.0) {
			// Otherwise, if the value of diag is not null, copy it in the
			// correct position to mantain the order
			j = b.dim;
			while ((j > 0) && (vpos[j-1] > (i+index3))) { 
				vpos[j] = vpos[j-1]; vval[j] = vval[j-1]; j--; 
			}
			vpos[j] = (i+index3); vval[j] = diag; b.dim++;
		}
		// Update the values to the next iteration
		vpos += b.dim; vval += b.dim; *(vlen++) = b.dim;
	}
	// Knowing the lengths, it is easy to obtain the vector vptr
	TransformLengthtoHeader (dst->vptr, dst->dim1);
	// Remove the permutations applied to make the addition
	reloj (&te1, &tu1);
	if (perm_cols) {
		PermuteColsSparseMatrix (src2, index2, iperm2); PermuteColsSparseMatrix (src1, index1, iperm1);
	}
	reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);
}

/*******************************************************************/
 
// This routine computes the number of elements which appears in both sparse matrices.
// The vector iperm1 contains the inverse permutation of perm1, in which the value
// -1 indicates that the correspondient position doesn't exist in perm1.
// The same for iperm2 and perm2.
// The parameters index1 and index2 indicate, respectivaly, if 0-indexing 
// or 1-indexing is used to store the different sparse matrices.
long int IntersectSparseMatrices (SparseMatrix src1, int index1, int *perm1, int *iperm1,
																	SparseMatrix src2, int index2, int *perm2, int *iperm2,
																	ptr_SparseMatrix dst, int index3) {
	int i;
	int *vlen = dst->vptr+1, *vpos = dst->vpos + *(dst->vptr) - index3;
	double diag, *vval = dst->vval + *(dst->vptr) - index3;
	SparseVector b1, b2, b;
	int    *pp1 = src1.vptr - index1, *pp2 = src2.vptr - index2;
	int    *pi1 = src1.vpos - index1, *pi2 = src2.vpos - index2;
	double *pd1 = src1.vval - index1, *pd2 = src2.vval - index2;
	long int Jcont = 0;

	// Apply the permutations required to make the addition
	PermuteColsSparseMatrix (src1, index1, perm1); PermuteColsSparseMatrix (src2, index2, perm2);

	for (i=0; i<dst->dim1; i++) {
		// Initialization of the elements on which the addition will be made
		diag = 0.0; b.dim = 0; b.vpos = vpos; b.vval = vval;
		if ((iperm1[i] != -1) && (iperm2[i] != -1)) {
			// If the row appears in both matrices, the addition has to be made
			// If the MSR is used, the diag is changed
			if (src1.vptr == src1.vpos) diag += pd1[iperm1[i]];
			if (src2.vptr == src2.vpos) diag += pd2[iperm2[i]];
			// The sparse vector b1 is made from the i-th row of src1
			b1.dim = pp1[iperm1[i]+1] - pp1[iperm1[i]];
			b1.vpos = pi1 + pp1[iperm1[i]]; b1.vval = pd1 + pp1[iperm1[i]];
			// The sparse vector b2 is made from the i-th row of src2
			b2.dim = pp2[iperm2[i]+1] - pp2[iperm2[i]];
			b2.vpos = pi2 + pp2[iperm2[i]]; b2.vval = pd2 + pp2[iperm2[i]];
			// The addition is computed
			Jcont += IntersectSparseVectors (b1, index1, b2, index2);
		} else if (iperm1[i] != -1) {
			// If the row only in matrix src1, it is copied
			// If the MSR is used, the diag is changed
			if (src1.vptr == src1.vpos) diag += pd1[iperm1[i]];
			// The sparse vector b1 is made from the i-th row of src1
			b1.dim = pp1[iperm1[i]+1] - pp1[iperm1[i]];
			b1.vpos = pi1 + pp1[iperm1[i]]; b1.vval = pd1 + pp1[iperm1[i]];
		} else { 
			// If the row only in matrix src2, it is copied
			// If the MSR is used, the diag is changed
			if (src2.vptr == src2.vpos) diag += pd2[iperm2[i]];
			// The sparse vector b2 is made from the i-th row of src2
			b2.dim = pp2[iperm2[i]+1] - pp2[iperm2[i]];
			b2.vpos = pi2 + pp2[iperm2[i]]; b2.vval = pd2 + pp2[iperm2[i]];
		}

		// Update the values to the next iteration
		vpos += b.dim; vval += b.dim; *(vlen++) = b.dim;
	}
	// Knowing the lengths, it is easy to obtain the vector vptr
	TransformLengthtoHeader (dst->vptr, dst->dim1);
	// Remove the permutations applied to make the addition
	PermuteColsSparseMatrix (src2, index2, iperm2); PermuteColsSparseMatrix (src1, index1, iperm1);

	return Jcont;
}

/*******************************************************************/

