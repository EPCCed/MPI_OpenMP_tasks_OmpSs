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
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectors.h"
#include "SparseMatrices.h"

void CreateSparseMatrix (ptr_SparseMatrix spr, int numR, int numC, int numE, int msr)
	{
		spr->dim1 = numR; spr->dim2 = numC; 
		CreateInts (&(spr->vptr), numE+numR+1);
		*(spr->vptr) = ((msr)? (numR+1): 0);
		spr->vpos = spr->vptr + ((msr)? 0: (numR+1));
		CreateDoubles (&(spr->vval), numE+(numR+1)*msr);
	}

void ReallocSparseMatrix (ptr_SparseMatrix spr, int msr) 
	{
		int numR = spr->dim1, numE = spr->vptr[numR];
	
		ReallocInts (&(spr->vptr), numE+numR+1);
		*(spr->vptr) = ((msr)? (numR+1): 0);
		spr->vpos = spr->vptr + ((msr)? 0: (numR+1));
		ReallocDoubles (&(spr->vval), numE+(numR+1)*msr);
	}

void ConvertSparseMatrix (SparseMatrix spr, int lngs)
	{
		int i;

		if ((spr.dim1 > 0) && (spr.dim2 > 0))
			if ((lngs & 2) != (lngs & 1)) {
				if (lngs & 1)
					{
						for (i=0; i<=spr.dim1; i++) spr.vptr[i]--;
						for (i=spr.vptr[0]; i<spr.vptr[spr.dim1]; i++) spr.vpos[i]--;
					}
				else
					{
						for (i=spr.vptr[0]; i<spr.vptr[spr.dim1]; i++) spr.vpos[i]++;
						for (i=0; i<=spr.dim1; i++) spr.vptr[i]++;
					}
			}
	}

void CopySparseMatrices (SparseMatrix src, ptr_SparseMatrix dst, int lngs)
	{
		int lng1 = (lngs & 1), lng2 = ((lngs & 2)>>1), off11 = *(src.vptr) , off21 = *(dst->vptr);
		int off12 = off11 - lng1 , off22 = off21 - lng2;
		int i, j, ns, nd;
		int *pps1, *pps2, *pis, *ppd1, *ppd2, *pid, *pi1;
		double *pds, *pdd, *pd1, *pd3;


		if (((src.vptr == src.vpos) && (dst->vptr == dst->vpos)) ||
				((src.vptr != src.vpos) && (dst->vptr != dst->vpos)))
			{
				if ((src.vptr == src.vpos) && (dst->vptr == dst->vpos))
					CopyDoubles (src.vval, dst->vval, src.dim1);
				CopyShiftInts (src.vptr+1, dst->vptr+1, src.dim1, off21-off11);
				CopyShiftInts (src.vpos + off12, dst->vpos+off22, src.vptr[src.dim1]-off11, lng2-lng1);
				CopyDoubles (src.vval + off12, dst->vval + off22, src.vptr[src.dim1] - off11);
			}
		else if (dst->vptr == dst->vpos) 
			{
				pps1 = src.vptr ; pps2 = pps1+1; ns = (*pps2-*pps1);
				ppd1 = dst->vptr; ppd2 = ppd1+1;
				pis =	src.vpos+*pps1; pds =	src.vval+*pps1;
				pid = dst->vpos+*ppd1; pdd = dst->vval+*ppd1;
 				pd3 = dst->vval;
				for (i=0; i<src.dim1; i++)
					{
						nd = ns; j = ns; pi1 = pis; pd1 = pds;
						while ((j > 0) && (*pi1 < i))
							{ pi1++; j--; }
						if ((j > 0) && (*pi1 == i))
							{
								CopyInts		(pis, pid, ns-j);
								CopyDoubles (pds, pdd, ns-j);
								pd1 += (ns-j); pid += (ns-j); pdd += (ns-j);
								*(pd3++) = *(pd1++); pi1++; j--;
								CopyInts		(pi1, pid, j);
								CopyDoubles (pd1, pdd, j);
								pid += j; pdd += j; nd--;
							}
						else
							{
								*(pd3++) = 0.0;
								CopyInts		(pis, pid, ns);
								CopyDoubles (pds, pdd, ns);
								pid += ns; pdd += ns;
							}
						pis += ns; pds += ns; 
						pps1 = (pps2++); ns = (*pps2-*pps1);
						*ppd2 = *ppd1 + nd; ppd1 = (ppd2++); 
					}
			}
		else  
			{
				pps1 = src.vptr ; pps2 = pps1+1; ns = (*pps2-*pps1);
				ppd1 = dst->vptr; ppd2 = ppd1+1;
				pis =	src.vpos+*pps1; pds =	src.vval+*pps1;
				pid = dst->vpos+*ppd1; pdd = dst->vval+*ppd1;
 				pd3 = src.vval;
				for (i=0; i<src.dim1; i++)
					{
						if (*pd3 != 0.0)
							{
								nd = ns; j = ns; pi1 = pis; pd1 = pds;
								while ((j > 0) && (*pi1 < i))
									{ pi1++; j--; }
								CopyInts		(pis, pid, ns-j);
								CopyDoubles (pds, pdd, ns-j);
								pd1 += (ns-j); pid += (ns-j); pdd += (ns-j);
								*(pid++) = i; *(pdd++) = *pd3;
								CopyInts		(pi1, pid, j);
								CopyDoubles (pd1, pdd, j);
								pid += j; pdd += j; nd++;
							}
						else
							{
								CopyInts		(pis, pid, ns);
								CopyDoubles (pds, pdd, ns);
								nd = ns; pid += ns; pdd += ns;
							}
						pd3++; pis += ns; pds += ns; 
						pps1 = (pps2++); ns = (*pps2-*pps1);
						*ppd2 = *ppd1 + nd; ppd1 = (ppd2++); 
					}
			}
	}



void PrintSparseMatrix (SparseMatrix spr, int CorF)
	{
		int i, j;

		if (spr.vptr == spr.vpos)
			{
				printf ("Diagonals : \n ");
				for (i=0; i<spr.dim1; i++) printf ("%f ", spr.vval[i]); printf ("\n");
			}

		printf ("Pointers: \n ");
		if (spr.dim1 > 0)
			for (i=0; i<=spr.dim1; i++) printf ("%d ", spr.vptr[i]); printf ("\n");

		printf ("Values: \n");
		for (i=0; i<spr.dim1; i++)
			{ printf (" Row %d --> ", i+CorF);
			for (j=(spr.vptr[i]-CorF); j<(spr.vptr[i+1]-CorF); j++)
				printf ("(%d,%f) ", spr.vpos[j], spr.vval[j]); 
			printf ("\n");	}
		printf ("\n");
	}

void FPrintSparseMatrix (FILE *f, SparseMatrix spr, int CorF)
	{
		int i, j;

		if (spr.vptr == spr.vpos)
			{
				fprintf (f, "Diagonals : \n ");
				for (i=0; i<spr.dim1; i++) fprintf (f, "%f ", spr.vval[i]); fprintf (f, "\n");
			}

		fprintf (f, "Pointers: \n ");
		if (spr.dim1 > 0)
			for (i=0; i<=spr.dim1; i++) fprintf (f, "%d ", spr.vptr[i]); fprintf (f, "\n");

		fprintf (f, "Values: \n");
		for (i=0; i<spr.dim1; i++)
			{ fprintf (f, " Row %d --> ", i+CorF);
			for (j=spr.vptr[i]-CorF; j<spr.vptr[i+1]-CorF; j++)
				fprintf (f, "(%d,%f) ", spr.vpos[j], spr.vval[j]); 
			fprintf (f, "\n");	}
		fprintf (f, "\n");
	}

void RemoveSparseMatrix (ptr_SparseMatrix spr)
	{
		spr->dim1 = -1; spr->dim2 = -1; 
		RemoveInts (&(spr->vptr));
		RemoveDoubles (&(spr->vval)); 
	}

void GetDiagonalSparseMatrix (SparseMatrix spr, double *diag)
	{
		int i, j, dim = (spr.dim1 < spr.dim2)? spr.dim1 : spr.dim2;
		int *pp1 = NULL, *pp2 = NULL, *pi1 = NULL; 
		double *pd1 = NULL, *pd2 = diag;

		if (spr.vptr == spr.vpos)
			CopyDoubles (spr.vval, diag, spr.dim1);
		else
			{
				pp1 = spr.vptr; pp2 = pp1+1; j = (*pp2-*pp1);
				pi1 = spr.vpos+*pp1; pd1 = spr.vval+*pp1; 
				for (i=0; i<dim; i++)
					{
						while ((j > 0) && (*pi1 < i))
							{ pi1++; pd1++; j--; }
						*(pd2++) = ((j > 0) && (*pi1 == i))? *pd1: 0.0;
						pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1);
					}
			}
	}

void GetDiagonalSparseMatrixDspls (SparseMatrix spr, double *diag, int dspls)
	{
		int i, j, dim = (spr.dim1 < spr.dim2)? spr.dim1 : spr.dim2;
		int *pp1 = NULL, *pp2 = NULL, *pi1 = NULL; 
		double *pd1 = NULL, *pd2 = diag;

		if (spr.vptr == spr.vpos)
			CopyDoubles (spr.vval+dspls, diag, spr.dim1-dspls);
		else
			{
				pp1 = spr.vptr+dspls; pp2 = pp1+1; j = (*pp2-*pp1);
				pi1 = spr.vpos+*pp1; pd1 = spr.vval+*pp1; 
				for (i=dspls; i<dim; i++)
					{
						while ((j > 0) && (*pi1 < i))
							{ pi1++; pd1++; j--; }
						*(pd2++) = ((j > 0) && (*pi1 == i))? *pd1: 0.0;
						pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1);
					}
			}
	}

void HeapSortSparseMatrix (SparseMatrix spr, SparseVector_Sort_func funcSort)
	{
		int i, j; 
		int *pp1 = spr.vptr, *pp2 = NULL, *pi1 = spr.vpos + *pp1;
		double *pd1 = spr.vval + *pp1;
		SparseVector baux;

		pp2 = pp1+1; j = (*pp2-*pp1);
		for (i=0; i<spr.dim1; i++)
			{
				baux.dim = j; baux.vpos = pi1; baux.vval = pd1;
				HeapSortSparseVector (baux, funcSort);
				pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1); 
			}
	}

void PermuteColsSparseMatrix (SparseMatrix spr, int *perm)
	{
		int i, j; 
		int *pp1 = spr.vptr, *pp2 = NULL, *pi1 = spr.vpos + *pp1;
		double *pd1 = spr.vval + *pp1;
		SparseVector baux;

		pp2 = pp1+spr.dim1; j = (*pp2-*pp1); 
		PermuteInts (pi1, perm, 0, j); 

		pp2 = pp1+1; j = (*pp2-*pp1);
		for (i=0; i<spr.dim1; i++)
			{
				baux.dim = j; baux.vpos = pi1; baux.vval = pd1;
				HeapSortSparseVector (baux, HeapSortSparseVector_IncrPos);
				pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1); 
			}
	}

void PermuteColsWithNegSparseMatrix (SparseMatrix spr, int *perm)
	{
		int i, j, k1 = 0, k2 = 0; 
		int *pp1 = spr.vptr, *pp2 = NULL, *pi1 = spr.vpos + *pp1, *pi2 = NULL;
		double *pd1 = spr.vval + *pp1;
		SparseVector baux;

		pp2 = pp1+spr.dim1; j = (*pp2-*pp1); if (perm != NULL) PermuteInts (pi1, perm, 0, j); 

		pp2 = pp1+1; j = (*pp2-*pp1);
		for (i=0; i<spr.dim1; i++)
			{
				baux.dim = j; baux.vpos = pi1; baux.vval = pd1;
				HeapSortSparseVector (baux, HeapSortSparseVector_IncrPosWithNeg);

				k1 = j; pi2 = pi1+j; while ((k1 > 0) && (*(--pi2) < 0)) k1--;
				if (k2 > 0)
					{
						CopyInts (pi1, pi1-k2, k1);
						CopyDoubles (pd1, pd1-k2, k1);
						*pp1 -= k2;
					}
				k2 += (j - k1); pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1);
			}
		*pp1 -= k2;
	}

void GetGraphWithDiagSparseMatrix (SparseMatrix spr, int *vptr, int *vpos)
	{
		int n = spr.dim1, i, j1, j2, dim, _bool, diag = 1;
		int *ptrSpr = NULL, *pi1 = spr.vptr, *pi2 = pi1 + 1, *pi3 = vpos;

		CreateInts (&ptrSpr, n); CopyInts (pi1, ptrSpr, n);
		*vptr = 0;
		for (i=0; i<n; i++)
			{
				j1 =	*pi1; j2 = 0; dim = 0; diag = 1;
				while ((j1 < *pi2) && (j2 < n))
					{
						_bool = 1;
						while ((j2 < n) && _bool)
							if (ptrSpr[j2] < spr.vptr[j2+1]) 
								if (spr.vpos[ptrSpr[j2]] < i) ptrSpr[j2]++;
								else if (spr.vpos[ptrSpr[j2]] > i) j2++;
								else _bool = 0; 
							else j2++; 
						if (_bool == 0)
							{
								while ((j1 < *pi2) && (spr.vpos[j1] < j2))
									{
										if (diag && (spr.vpos[j1] > i)) 
											{ *(pi3++) = i; dim++; diag = 0; }
										else diag &= (i != spr.vpos[j1]);
										*(pi3++) = (spr.vpos[j1++]); dim++; 
									}
								if (j1 < *pi2)
									{
										if (spr.vpos[j1] == j2)
											{
												if (diag && (j2 > i)) 
													{ *(pi3++) = i; dim++; diag = 0; }
												else diag &= (i != spr.vpos[j1]);
												*(pi3++) = (spr.vpos[j1++]); dim++; ptrSpr[j2]++; j2++; 
											}
										else 
											{
												if (diag && (j2 > i)) 
													{ *(pi3++) = i; dim++; diag = 0; }
												*(pi3++) = j2; ptrSpr[j2]++;	dim++; 
											}
									}
								else 
									{
										if (diag && (j2 > i)) 
											{ *(pi3++) = i; dim++; diag = 0; }
										else diag &= (i != j2);
										*(pi3++) = j2; ptrSpr[j2]++;	dim++; 
									}
							}
					}
				while (j1 < *pi2) 
					{
						if (diag && (spr.vpos[j1] > i)) 
							{ *(pi3++) = i; dim++; diag = 0; }
						else diag &= (i != spr.vpos[j1]);
						*(pi3++) = (spr.vpos[j1++]); dim++; 
					}
				while (j2 < n)
					if (ptrSpr[j2] < spr.vptr[j2+1]) 
						if (spr.vpos[ptrSpr[j2]] < i) ptrSpr[j2]++;
						else if (spr.vpos[ptrSpr[j2]] > i) j2++;
						else 
							{
								if (diag && (j2 > i)) 
									{ *(pi3++) = i; dim++; diag = 0; }
								else diag &= (i != j2);
								*(pi3++) = j2; ptrSpr[j2]++;	dim++; 
							}
					else j2++;
				if (diag) { *(pi3++) = i; dim++; diag = 0; }
				pi1 = pi2; pi2++; *(vptr+1) = *vptr + dim; vptr++;
			}

		RemoveInts (&ptrSpr);
	}

void GetGraphSparseMatrix (SparseMatrix spr, int *vptr, int *vpos)
	{
		int n = spr.dim1, i, j1, j2, dim, _bool;
		int *ptrSpr = NULL, *pi1 = spr.vptr, *pi2 = pi1 + 1, *pi3 = vpos;
		FILE *f;

		CreateInts (&ptrSpr, n); CopyInts (pi1, ptrSpr, n);
		*vptr = 0;
		for (i=0; i<n; i++)
			{
        f = fopen ("Trace/METIS", "a"); fprintf (f, "%d \n", i); fclose (f);
				j1 =	*pi1; j2 = 0; dim = 0; 
				while ((j1 < *pi2) && (j2 < n))
					{
						_bool = 1;
						while ((j2 < n) && _bool)
							if (ptrSpr[j2] < spr.vptr[j2+1]) 
								if (spr.vpos[ptrSpr[j2]] < i) ptrSpr[j2]++;
								else if (spr.vpos[ptrSpr[j2]] > i) j2++;
								else _bool = 0; 
							else j2++; 
						if (_bool == 0)
							{
								while ((j1 < *pi2) && (spr.vpos[j1] < j2))
									if (spr.vpos[j1] != i) 
										{ *(pi3++) = (spr.vpos[j1++]); dim++; }
									else j1++;
								if (j1 < *pi2)
									{
										if (spr.vpos[j1] == j2)
											{
												if (spr.vpos[j1] != i) 
													{ *(pi3++) = (spr.vpos[j1++]); dim++; ptrSpr[j2]++; j2++; }
												else { j1++; ptrSpr[j2]++; j2++; }
											}
										else 
											{
												if (j2 != i) 
													{ *(pi3++) = j2; ptrSpr[j2]++;	dim++; }
												else { ptrSpr[j2]++; j2++; }
											}
									}
								else 
									{
										if (j2 != i) 
											{ *(pi3++) = j2; ptrSpr[j2]++;	dim++; }
										else { ptrSpr[j2]++; j2++; }
									}
							}
					}
				while (j1 < *pi2) 
					{
						if (spr.vpos[j1] != i)
							{ *(pi3++) = (spr.vpos[j1++]); dim++; }
						else j1++;
					}
				while (j2 < n)
					if (ptrSpr[j2] < spr.vptr[j2+1]) 
						if (spr.vpos[ptrSpr[j2]] < i) ptrSpr[j2]++;
						else if (spr.vpos[ptrSpr[j2]] > i) j2++;
						else 
							{
								if (j2 != i) 
									{ *(pi3++) = j2; ptrSpr[j2]++;	dim++; }
								else { ptrSpr[j2]++; j2++; }
							}
					else j2++;
				pi1 = pi2; pi2++; *(vptr+1) = *vptr + dim; vptr++;
			}

		RemoveInts (&ptrSpr);
	}


void GetGraphSparseMatrix2 (SparseMatrix spr, int *vptr, int *vpos)
	{
		int n = spr.dim1, row_id, count;
		int *pi1 = spr.vptr, *pi2 = pi1 + 1, *pi3 = vptr + n, *pi4 = vptr + 1, *pi5 ;
		int * tmp_pointers = NULL;

		CreateInts (&tmp_pointers, n+1);

		//Calculate each spr row size on vptr[1:n]
		*(pi4 - 1) = 0;
		while (pi4 <= pi3)
			{
				*pi4 = (*pi2 - *pi1) - 1;
				pi1 = pi2++;
				pi4++;
			}


		pi1 = spr.vptr;
		pi2 = pi1 + 1;
		pi3 = spr.vptr + n ;

		while (pi2 <= pi3)
			{
				//Skip diagonal entry
				pi4 = spr.vpos + *pi1 + 1;
				pi5 = spr.vpos + *pi2;
	
				while (pi4 < pi5)
					{
						(*(vptr + (*pi4 + 1)))++;
						pi4++;
					}
				pi1 = pi2++;
			}

		TransformLengthtoHeader(vptr,n);
		CopyInts (vptr, tmp_pointers, n+1);

		pi1 = spr.vptr ;
		pi2 = pi1 + 1 ;
		pi3 = spr.vptr + n ;

		while (pi2 <= pi3)
			{
				row_id = pi1 - spr.vptr ;
				count = ( (*pi2 - *pi1) -1 ) ;
		
				//Skip diagonal entry
				pi4 = spr.vpos + *pi1 + 1 ;

				pi5 = spr.vpos + *pi2 ;

				//Copy adjacent nodes present in the original graph
				CopyInts (pi4, vpos+tmp_pointers[row_id], count);
				tmp_pointers[row_id] += count;

				//Copy symmetrized adjacent nodes
				while (pi4 < pi5)
					{
						*(vpos + tmp_pointers[*pi4]) = row_id;
						tmp_pointers[*pi4]++;
						pi4++;
					}

					pi1 = pi2++;
			}
		RemoveInts(&tmp_pointers);

	}

void GetGraphSparseMatrix3 (SparseMatrix spr, int indexS, int *vptr, int *vpos, int indexG) {
	int n = spr.dim1;
	int * tmp_pointers = NULL;
	int *pp1 = NULL, *pp2 = NULL, *pp3 = NULL, *pp4 = NULL, *pp5 = NULL;
	int i, j, dim, indexGS = indexG - indexS;

	CreateInts (&tmp_pointers, n+1);

	//Calculate each spr row size on vptr[1:n]
	*vptr = indexG; InitInts (vptr+1, n, 0, 0);

	pp1 = spr.vptr; pp3 = spr.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = vptr+1-indexS;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1); pp4[i] += dim;
		for (j=0; j<dim; j++) {
			if (*pp3 == i) pp4[i]--; 
			else        pp4[*pp3]++;
			pp3++;
		}
		pp1 = pp2++; 
	}

	TransformLengthtoHeader (vptr, n);
	CopyInts (vptr, tmp_pointers, n+1);

	pp1 = spr.vptr; pp3 = spr.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = vpos-indexG; pp5 = tmp_pointers - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1);
		for (j=0; j<dim; j++) {
			if (*pp3 != i) {
				pp4[pp5[i]++] = *pp3+indexGS; 
				pp4[pp5[*pp3]++] = i+indexGS;
			}
			pp3++;
		}
		pp1 = pp2++;
	}

	RemoveInts(&tmp_pointers);
}

void ProdSparseMatrixVector (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos+*pp1;
		double *pd1 = spr.vval+*pp1, *pd2 = res;
		SparseVector baux;

		if (spr.vptr == spr.vpos)
			VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);

		j = (*pp2-*pp1);
		for (i=0; i<spr.dim1; i++)
			{
				baux.dim = j; baux.vpos = pi1; baux.vval = pd1;
				*(pd2++) += DotSparseVector (baux, vec);
				pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1); 
			}
	}

void ProdSparseMatrixVector2 (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		double aux;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++)
					aux += spr.vval[j] * vec[spr.vpos[j]];
				res[i] += aux;
			}
	}

void ProdSparseMatrixVector3 (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		double aux;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++)
					aux += spr.vval[j] * vec[spr.vpos[j]];
				res[i] = aux;
			}
	}

void AddSparseMatrices (SparseMatrix src1, int *perm1, int *iperm1,
												SparseMatrix src2, int *perm2, int *iperm2,
												ptr_SparseMatrix dst)
	{
		int i, j, k;
		int *vlen = dst->vptr+1, *vpos = dst->vpos + *(dst->vptr);
		double diag, *vval = dst->vval + *(dst->vptr);
		SparseVector b1, b2, b;

		PermuteColsSparseMatrix (src1, perm1); PermuteColsSparseMatrix (src2, perm2);

		for (i=0; i<dst->dim1; i++)
			{
				diag = 0.0;
				b.dim = 0; b.vpos = vpos; b.vval = vval;
				if ((iperm1[i] != -1) && (iperm2[i] != -1))
					{
						if (src1.vptr == src1.vpos) diag += src1.vval[iperm1[i]];
						if (src2.vptr == src2.vpos) diag += src2.vval[iperm2[i]];
						b1.dim = src1.vptr[iperm1[i]+1] - src1.vptr[iperm1[i]];
						b1.vpos = src1.vpos + src1.vptr[iperm1[i]];
						b1.vval = src1.vval + src1.vptr[iperm1[i]];
						b2.dim = src2.vptr[iperm2[i]+1] - src2.vptr[iperm2[i]];
						b2.vpos = src2.vpos + src2.vptr[iperm2[i]];
						b2.vval = src2.vval + src2.vptr[iperm2[i]];
						AddSparseVectors (b1, b2, &b);
					}
				else if (iperm1[i] != -1) 
					{
						if (src1.vptr == src1.vpos) diag += src1.vval[iperm1[i]];
						b1.dim = src1.vptr[iperm1[i]+1] - src1.vptr[iperm1[i]];
						b1.vpos = src1.vpos + src1.vptr[iperm1[i]];
						b1.vval = src1.vval + src1.vptr[iperm1[i]];
						CopySparseVectors (b1, &b);
					}
				else if (iperm2[i] != -1) 
					{
						if (src2.vptr == src2.vpos) diag += src2.vval[iperm2[i]];
						b2.dim = src2.vptr[iperm2[i]+1] - src2.vptr[iperm2[i]];
						b2.vpos = src2.vpos + src2.vptr[iperm2[i]];
						b2.vval = src2.vval + src2.vptr[iperm2[i]];
						CopySparseVectors (b2, &b);
					}

				if ((src1.vptr != src1.vpos) || (src2.vptr != src2.vpos))
					{
						k = FindSortInt (i, vpos, b.dim, HeapSortInts_IncrPos); 
						if (k != -1) 
							{
								if (dst->vptr == dst->vpos)	
									{
										*(dst->vval+i) = diag + *(vval+k); b.dim--;
										for (j=k; j<b.dim; j++) 
											{ vpos[j] = vpos[j+1]; vval[j] = vval[j+1]; }
									}
								else if (diag != 0.0)
									*(vval+k) += diag;
							}
						else if (dst->vptr == dst->vpos) 
							*(dst->vval+i) = diag;
						else if (diag != 0.0)
							{
								j = b.dim;
								while ((j > 0) && (vpos[j-1] > i))
									{ vpos[j] = vpos[j-1]; vval[j] = vval[j-1]; j--; }
								vpos[j] = i; vval[j] = diag; b.dim++;
							}
					}
				else if (dst->vptr == dst->vpos)	
					*(dst->vval+i) = diag;
				else if (diag != 0.0)
					{
						j = b.dim;
						while ((j > 0) && (vpos[j-1] > i))
							{ vpos[j] = vpos[j-1]; vval[j] = vval[j-1]; j--; }
						vpos[j] = i; vval[j] = diag; b.dim++;
					}

				vpos += b.dim; vval += b.dim; *(vlen++) = b.dim;
			}
		TransformLengthtoHeader (dst->vptr, dst->dim1);
		PermuteColsSparseMatrix (src2, iperm2); PermuteColsSparseMatrix (src1, iperm1);
	}
