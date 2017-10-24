#include <stdio.h>
#include <stdlib.h>
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectors.h"
#include "SparseMatrices.h"

void DesymmetrizeSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst) {
	int i, j, row, col, pos1, pos2;
	int n = src.dim1, nnz = 0;
	int *sizes = NULL;

	CreateInts (&sizes, n);
	InitInts (sizes, n, 0, 0);
	for (i=0; i<n; i++) {
		for (j=src.vptr[i]; j<src.vptr[i+1]; j++) {
			sizes[i]++; nnz++;
			if (src.vpos[j] != i) {
				sizes[src.vpos[j]]++; nnz++;
			}
		}
	}
	
	CreateSparseMatrix (dst, n, n, nnz, 0);
	CopyInts (sizes, (dst->vptr)+1, n);
	dst->vptr[0] = 0; TransformLengthtoHeader (dst->vptr, n);
	CopyInts (dst->vptr, sizes, n);
	for (i=0; i<n; i++) {
		for (j=src.vptr[i]; j<src.vptr[i+1]; j++) {
			row = i; pos1 = sizes[row]; 
			dst->vpos[pos1] = src.vpos[j];
			dst->vval[pos1] = src.vval[j];
			sizes[row]++;
			if (src.vpos[j] != i) {
				col = src.vpos[j]; pos2 = sizes[col];
				dst->vpos[pos2] = row;
				dst->vval[pos2] = src.vval[j];
				sizes[col]++;
			}
		}
	}
	RemoveInts (&sizes);
}

void RestrictUpperSparseMatrices (SparseMatrix src, ptr_SparseMatrix dst) {
	int i, j, k;
	int n = src.dim1, nnz = src.vptr[n];
	int *sizes = NULL;

	HeapSortSparseMatrix (src, HeapSortSparseVector_IncrPos);
	CreateInts (&sizes, n+1); CopyInts (src.vptr, sizes, n+1);

	CreateSparseMatrix (dst, n, n, nnz, 0);
	k = 0;
	for (i=0; i<n; i++) {
		dst->vptr[i] = k; j = i+1;
		while ((j < n) && (sizes[i] < src.vptr[i+1])) {
			while ((j < n) && (src.vpos[sizes[j]] != i)) j++;
			if (j < n) {
				while ((sizes[i] < src.vptr[i+1]) && (j > src.vpos[sizes[i]])) {
					dst->vpos[k] =  src.vpos[sizes[i]];
					dst->vval[k] =  src.vval[sizes[i]];
					k++; sizes[i]++;
				}
				if (j == src.vpos[sizes[i]]) {
					dst->vpos[k] =  src.vpos[sizes[i]] + src.vpos[sizes[j]];
					dst->vval[k] =  src.vval[sizes[i]] + src.vval[sizes[j]];
					k++; sizes[i]++; sizes[j]++; j++;
				} else {
					dst->vpos[k] =  j;
					dst->vval[k] =  src.vval[sizes[j]];
					k++; sizes[j]++; j++;
				}
			}
		}
		while (j < n) {
			if (src.vpos[sizes[j]] == i) {
				dst->vpos[k] =  j;
				dst->vval[k] =  src.vval[sizes[j]];
				k++; sizes[j]++; 
			}
			j++;
		}
		while (sizes[i] < src.vptr[i+1]) {
			dst->vpos[k] =  src.vpos[sizes[i]];
			dst->vval[k] =  src.vval[sizes[i]];
			k++; sizes[i]++;
		}
	}
	dst->vptr[n] = k;
	
	RemoveInts (&sizes);
}

void ProdSymSparseMatrixVector (SparseMatrix spr, double *vec, double *res)
	{
		int i, j, k;
		int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos+*pp1;
		double *pd1 = spr.vval+*pp1, *pd2 = vec, *pd3 = res;;
		SparseVector baux;

		if (spr.vptr == spr.vpos)
			{
				VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);
				j = (*pp2-*pp1);
				for (i=0; i<spr.dim1; i++)
					{
						baux.dim = j; baux.vpos = pi1; baux.vval = pd1;
						AxpySparseVector (baux, *(pd2++), res);
						*(pd3++) += DotSparseVector (baux, vec);
						pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1);
					}
			}
		else
			{
				j = (*pp2-*pp1);
				for (i=0; i<spr.dim1; i++)
					{
						baux.dim = j; baux.vpos = pi1; baux.vval = pd1;
						k = FindSortInt (i, pi1, j, HeapSortInts_IncrPos); // printf ("KKKKK = %d\n", k);
						if (k != -1) *pd3 -= *(pd1+k) * *pd2;
						AxpySparseVector (baux, *(pd2++), res);
						*(pd3++) += DotSparseVector (baux, vec);
						pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1); 
					}
			}
	}

void ProdSymSparseMatrixVector2 (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		double aux;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++) {
					aux += spr.vval[j] * vec[spr.vpos[j]];
					if (spr.vpos[j] != i)
						res[spr.vpos[j]] += spr.vval[j] * vec[i];
				}
				res[i] += aux;
			}
	}

void ProdSymSparseMatrixVector3 (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		double aux;

		InitDoubles (res, spr.dim1, 0.0, 0.0);
		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++) {
					aux += spr.vval[j] * vec[spr.vpos[j]];
					if (spr.vpos[j] != i)
						res[spr.vpos[j]] += spr.vval[j] * vec[i];
				}
				res[i] += aux;
			}
	}

