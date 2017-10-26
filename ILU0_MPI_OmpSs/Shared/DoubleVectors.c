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
#include <string.h>
#include "DoubleVectors.h"

void CreateDoubles (double **vdbl, int dim)
	{
		if ((*vdbl = (double *) malloc (sizeof(double)*dim)) == NULL)
			{ printf ("Memory Error (CreateDoubles(%d))\n", dim); exit (1); }
	}

void ReallocDoubles (double **vdbl, int dim)
	{
		if ((*vdbl = (double *) realloc (*vdbl, sizeof(double)*dim)) == NULL)
			{ printf ("Memory Error (ReallocDoubles(%d))\n", dim); exit (1); }
	}

void InitDoubles (double *vdbl, int dim, double frst, double incr) 
	{
		int i; 
		double *pd = vdbl, num = frst;

		for (i=0; i<dim; i++) 
			{ *(pd++) = num; num += incr; }
	}
		
void InitRandDoubles (double *vdbl, int dim, double frst, double last) 
	{
		int i; 
		double *pd = vdbl, size = last - frst;

		for (i=0; i<dim; i++) 
			{ *(pd++) = frst + (size * (rand() / (RAND_MAX + 1.0))); }
	}
		
void CopyDoubles (double *src, double *dst, int dim)
	{ 
		memmove (dst, src, sizeof(double) * dim);
	}

void ScaleDoubles (double *vdbl, double scal, int dim) 
	{
		int i; 
		double *pd = vdbl;

		for (i=0; i<dim; i++) 
			*(pd++) *= scal;
	}

void AxpyDoubles (double alfa, double *vdbl1, double *vdbl2, int dim) 
	{
		int i; 
		double *pd1 = vdbl1, *pd2 = vdbl2;

		for (i=0; i<dim; i++) 
			*(pd2++) += (*(pd1++) * alfa);
	}

void XpayDoubles (double *vdbl1, double alfa, double *vdbl2, int dim) 
	{
		int i; 
		double *pd1 = vdbl1, *pd2 = vdbl2;

		for (i=0; i<dim; i++) 
			{ *pd2 = ((*pd2) * alfa) + (*(pd1++)); pd2++; }
	}

double DotDoubles (double *vdbl1, double *vdbl2, int dim) 
	{
		int i; 
		double *pd1 = vdbl1, *pd2 = vdbl2, res = 0.0;

		for (i=0; i<dim; i++) 
			res += (*(pd2++)) * (*(pd1++));

		return res;
	}

void GetDoubleFromString (char *string, double *pdbl, int dimC)
	{
		int j, k, exp, neg;
		double num, frac;
		char *pchar = string;

		j = 0; exp = 0; neg = 0; num = 0.0; frac = 1.0;
		while ((j < dimC) && ((*pchar < '0') || (*pchar > '9')) && 
					 (*pchar != '+') && (*pchar != '-') && (*pchar != '.')) { j++; pchar++; }
		if (j < dimC)
			{
				if ((*pchar == '+') || (*pchar == '-'))
					{ neg = (*pchar == '-'); j++; pchar++; }
				if (j < dimC)
					{
						if (*pchar != '.')
							while ((j < dimC) && (*pchar >= '0') && (*pchar <= '9'))
								{ num = num * 10 + (*pchar - 48); j++; pchar++; }
						if (j < dimC)
							{
								if (*pchar == '.')
									{
										j++; pchar++; 
										while ((j < dimC) && (*pchar >= '0') && (*pchar <= '9'))
											{ frac /= 10; num += (*pchar-48) * frac; j++; pchar++; }
									}
								if (neg) num = -num;
								if (j < dimC)
									{
										if ((*pchar == 'e') || (*pchar == 'E') || (*pchar == 'd') || (*pchar == 'D'))
											{
												neg = 0; j++; pchar++; 
												if (j < dimC)
													{
														if ((*pchar == '+') || (*pchar == '-'))
															{ neg = (*pchar == '-'); j++; pchar++; }
														if (j < dimC)
															{
																while ((j < dimC) && (*pchar >= '0') && 
																		 (*pchar <= '9'))
																	{ exp = exp*10 + (*pchar-48); j++; pchar++; }
																if (neg) exp = -exp;
																for (k=0; k<exp; k++) num *= 10;
																for (k=0; k>exp; k--) num /= 10;
															}
													}
											}
									}
							}
						else
 							if (neg) num = -num;
					}
			}
		*pdbl = num; 
	}

void GetDoublesFromString (char *string, double *vdbl, int dimN, int dimC)
	{
		int i;
		double *paux = vdbl;
		char *pchar = string;

		for (i=0; i<dimN; i++)
			{ GetDoubleFromString (pchar, (paux++), dimC); pchar += dimC; }
	}

int HeapSortPermDoubles_IncrPos (double n1, double n2, int p1, int p2)
	{ return ((n1 < n2) || ((n1 == n2) && (p1 < p2))); }

int HeapSortPermDoubles_DecrPos (double n1, double n2, int p1, int p2)
	{ return ((n1 > n2) || ((n1 == n2) && (p1 > p2))); }

int HeapSortPermDoubles_IncrVal (double n1, double n2, int p1, int p2)
	{ return (p1 < p2); }

int HeapSortPermDoubles_DecrVal (double n1, double n2, int p1, int p2)
	{ return (p1 > p2); }

void HeapSortPermDoubles_Swap (double *n1, double *n2, int *p1, int *p2)
	{
		double n;
		int p;
		n = *n1; *n1 = *n2; *n2 = n;
		p = *p1; *p1 = *p2; *p2 = p;
	}

void HeapSortPermDoubles_Adjust (double *vdbl, int *prm, int inic, int final, Doubles_SortPerm_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vdbl[j], vdbl[j+1], prm[j], prm[j+1]));
				_bool = funcSort (vdbl[k], vdbl[j], prm[k], prm[j]);
				if (_bool)
					{ 
						HeapSortPermDoubles_Swap (vdbl+k, vdbl+j, prm+k, prm+j); 
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortPermDoubles (double *vdbl, int *prm, int dim, Doubles_SortPerm_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortPermDoubles_Adjust (vdbl, prm, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortPermDoubles_Swap (vdbl, vdbl+(ind+1), prm, prm+(ind+1));
				HeapSortPermDoubles_Adjust(vdbl, prm, 0, ind, funcSort);
			}
	}

int HeapSortDoubles_IncrPos (double n1, double n2)
	{ return (n1 < n2); }

int HeapSortDoubles_DecrPos (double n1, double n2)
	{ return (n1 > n2); }

void HeapSortDoubles_Swap (double *n1, double *n2)
	{
		double n;
		n = *n1; *n1 = *n2; *n2 = n;
	}

void HeapSortDoubles_Adjust (double *vdbl, int inic, int final, Doubles_Sort_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vdbl[j], vdbl[j+1]));
				_bool = funcSort (vdbl[k], vdbl[j]);
				if (_bool)
					{ 
						HeapSortDoubles_Swap (vdbl+k, vdbl+j);
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortDoubles (double *vdbl, int dim, Doubles_Sort_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortDoubles_Adjust (vdbl, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortDoubles_Swap (vdbl, vdbl+(ind+1));
				HeapSortDoubles_Adjust(vdbl, 0, ind, funcSort);
			}
	}

int FindSortDouble (int num, double *vdbl, int dim, Doubles_Sort_func funcSort)
	{
		int _bool = 1, izq = 0, der = dim - 1, cnt = 0; 

		while ((_bool) && (izq <= der))
			{
				cnt = (izq + der) / 2;
				if (num == vdbl[cnt])
					_bool = 0;
				else if (funcSort (num, vdbl[cnt])) der = cnt - 1;
				else izq = cnt + 1;
			}
		if (_bool) cnt = -1;
		else if (cnt > 0)
			{
				if (funcSort (0,1))
					{ while ((cnt >= 0) && (vdbl[cnt] == num)) cnt--; cnt++; }
				else
					{ while ((cnt < dim) && (vdbl[cnt] == num)) cnt++; cnt--; }
			}
		return cnt;
	}
	
void VdivDoubles (double alfa, double *src, double *dst, int dim)
  {
    int i;

    for (i=0; i<dim; i++)
      { *dst = (alfa * *dst / *(src++)); dst++; }
  }

void VvecDoubles (double alfa, double *src1, double *src2, double beta, double *dst, int dim)
  {
    int i;

    for (i=0; i<dim; i++)
      { *dst = (beta * *dst) + (alfa * *(src1++) * *(src2++)); dst++; }
  }

double AddDoubles (double *vdbl, int dim) 
	{
		int i;
		double *pd = vdbl, aux = 0;

		for (i=0; i<dim; i++) { aux += *pd; pd++; }
		return aux;
	}

double AddPermuteDoubles (double *vdbl, int *perm, int index, int dim) 
	{
		int i, *pi = perm;
		double aux = 0, *pd = vdbl - index;
		
		for (i=0; i<dim; i++) { aux += pd[*pi]; pi++; }
		return aux;
	}

double MaxDoubles (double *vdbl, int dim) 
	{
		int i;
		double *pd = vdbl, aux = *pd;

		for (i=1; i<dim; i++) { if (*pd > aux) aux = *pd; pd++; }
		return aux;
	}

void CopyPermuteDoubles (double *src, int *perm, int index, double *dst, int dim)
	{
		int i, *pp = perm;
		double *pd1 = src - index, *pd2 = dst;
		
		for (i=0; i<dim; i++) { *pd2 = pd1[*pp]; pp++; pd2++; }
	}

void CopyInvPermuteDoubles (double *src, double *dst, int *perm, int index, int dim)
	{
		int i, *pp = perm;
		double *pd1 = src, *pd2 = dst - index;
		
		for (i=0; i<dim; i++) { pd2[*pp] = *pd1; pp++; pd1++; }
	}

void PrintDoubles (double *vdbl, int dim)
	{
		int i;
		double *pd = vdbl;

		for (i=0; i<dim; i++) printf ("%e ", *(pd++));
		printf ("\n");
	}

void FPrintDoubles (FILE *f, double *vdbl, int dim)
	{
		int i;
		double *pd = vdbl;

		for (i=0; i<dim; i++) fprintf (f, "%20.14e ", *(pd++));
		fprintf (f, "\n");
	}

void PrintFDoubles (double *vdbl, int dim, int f1, int f2) {
	char formato[10];
	int i;
	double *pd = vdbl;

	if (f1 <= 0)
		sprintf (formato, "%%e ");
	else if (f2 <= 0)
		sprintf (formato, "%%%de ", f1);
	else
		sprintf (formato, "%%%d.%de ", f1,f2);

	for (i=0; i<dim; i++) printf (formato, *(pd++));
	printf ("\n");
}

void FPrintFDoubles (FILE *f, double *vdbl, int dim, int f1, int f2) {
	char formato[10];
	int i, j, num;
	double *pd = vdbl;

	if (f1 <= 0)
		sprintf (formato, "%%e ");
	else if (f2 <= 0)
		sprintf (formato, "%%%de ", f1);
	else
		sprintf (formato, "%%%d.%de ", f1,f2);

	num = (f1 <= 0)? 10: (80 / f1); j = num; 
	for (i=0; i<dim; i++) {
		fprintf (f, formato, *(pd++));
		if ((--j) == 0) {
			fprintf (f, "\n"); j = num;
		}
	}		
	if (j < num) fprintf (f, "\n");
}

int ReadDoubles (char *file, double **vdbl) {
	int n = -1;
	double num;
	FILE *f = NULL;

	f = fopen (file, "r");
	if (f != NULL) {
		n = 0;
		while (fscanf (f, "%lf", &num) != EOF) n++;
		fclose (f);

		CreateDoubles (vdbl, n); n = 0;
		f = fopen (file, "r"); 
		while (fscanf (f, "%lf", &num) != EOF) 
			(*vdbl)[n++] = num;
		fclose (f);
	}

	return n;
}

void WriteDoubles (char *filename, double *vdbl, int dim) {
	FILE *f = NULL;

	printf ("Print %s of size %d\n", filename, dim);
	f = fopen (filename, "w");
	FPrintDoubles (f, vdbl, dim);
	fclose (f);
}

void WriteFDoubles (char *filename, double *vdbl, int dim, int f1, int f2) {
	FILE *f = NULL;

	printf ("Print %s of size %d\n", filename, dim);
	f = fopen (filename, "w");
	FPrintFDoubles (f, vdbl, dim, f1, f2);
	fclose (f);
}

void RemoveDoubles (double **vdbl)
	{ if (*vdbl != NULL) free (*vdbl); *vdbl = NULL; }

/****************************************************************************/

void CreateVectorPtrDoubles (matDoubles *mdbl, int numR)
	{
		if ((*mdbl = (matDoubles) malloc (sizeof(double *) * numR)) == NULL)
			{ printf ("Memory Error (CreateMatrixDoubles(%d))\n", numR); exit (1); }
	}

void CopyVectorPtrDoubles (matDoubles msrc, matDoubles mdst, int dim)
	{ memmove (mdst, msrc, sizeof(double *) * dim); }

void DesplVectorPtrDoubles (matDoubles mdbl, int dim, int dspl)
	{ 
		int i; 
		for (i=0; i<dim; i++) mdbl[i] += dspl; 
	}

void RemoveVectorPtrDoubles (matDoubles *mdbl)
	{
		if (*mdbl != NULL) free (*mdbl); 
		*mdbl = NULL;
	}

void CreateMatrixDoubles (matDoubles *mdbl, int numR, int numC)
	{
		int i;

		CreateVectorPtrDoubles (mdbl, numR);
		CreateDoubles (*mdbl, numR*numC);
		for (i=1; i<numR; i++) (*mdbl)[i] = (*mdbl)[i-1] + numC;
	}

void BuildMatrixDoubles (matDoubles mdbl, double *vec, int numR, int numC)
	{
		int i;

		*mdbl = vec;
		for (i=1; i<numR; i++) (*mdbl)[i] = (*mdbl)[i-1] + numC;
	}

void RemoveMatrixDoubles (matDoubles *mdbl)
	{
		if (*mdbl != NULL) RemoveDoubles (*mdbl); 
		RemoveVectorPtrDoubles (mdbl);
	}

