#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ScalarVectors.h"
#include "SparseVectors.h"

void CreateSparseVector (ptr_SparseVector b, int nz)
	{
		b->dim = 0;
		CreateInts (&(b->vpos), nz);
		CreateDoubles (&(b->vval), nz);
	}

void CopySparseVectors (SparseVector src, ptr_SparseVector dst)
	{ 
		dst->dim = src.dim;
		CopyInts (src.vpos, dst->vpos, src.dim);
		CopyDoubles (src.vval, dst->vval, src.dim);
	}

void PrintSparseVectors (SparseVector b)
	{
		int i, *pi = b.vpos;
		double *pd = b.vval;

		for (i=0; i<b.dim; i++)
			printf ("(%d,%f) ", *(pi++), *(pd++));
		printf ("\n");
	}

void FPrintSparseVectors (FILE *f, SparseVector b)
	{
		int i, *pi = b.vpos;
		double *pd = b.vval;

		for (i=0; i<b.dim; i++)
			fprintf (f, "(%d,%f) ", *(pi++), *(pd++));
		fprintf (f, "\n");
	}

void RemoveSparseVector (ptr_SparseVector b)
	{ b->dim = -1; RemoveDoubles (&(b->vval)); RemoveInts (&(b->vpos)); }

int HeapSortSparseVector_IncrPosWithNeg (int n1, int n2, double x1, double x2)
	{ return ((n2 == -1) || ((n1 != -1) && (n1 < n2))); }

int HeapSortSparseVector_DecrPosWithNeg (int n1, int n2, double x1, double x2)
	{ return ((n1 == -1) || ((n2 != -1) && (n1 > n2))); }

int HeapSortSparseVector_IncrPos (int n1, int n2, double x1, double x2)
	{ return (n1 < n2); }

int HeapSortSparseVector_DecrPos (int n1, int n2, double x1, double x2)
	{ return (n1 > n2); }

int HeapSortSparseVector_IncrVal (int n1, int n2, double x1, double x2)
	{ return (x1 < x2); }

int HeapSortSparseVector_DecrVal (int n1, int n2, double x1, double x2)
	{ return (x1 > x2); }

int HeapSortSparseVector_IncrAbsVal (int n1, int n2, double x1, double x2)
	{ return (fabs(x1) < fabs(x2)); }

int HeapSortSparseVector_DecrAbsVal (int n1, int n2, double x1, double x2)
	{ return (fabs(x1) > fabs(x2)); }

void HeapSortSparseVector_Swap (int *n1, int *n2, double *x1, double *x2)
	{
		int n; double x;
		n = *n1; *n1 = *n2; *n2 = n;
		x = *x1; *x1 = *x2; *x2 = x;
	}

void HeapSortSparseVector_Adjust (int *vpos, double *vval, int inic, int final, SparseVector_Sort_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vpos[j], vpos[j+1], vval[j], vval[j+1]));
				_bool = funcSort (vpos[k], vpos[j], vval[k], vval[j]);
				if (_bool)
					{ 
						HeapSortSparseVector_Swap (vpos+k, vpos+j, vval+k, vval+j); 
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortSparseVector (SparseVector b, SparseVector_Sort_func funcSort)
	{
		int i, ind, medio = b.dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortSparseVector_Adjust (b.vpos, b.vval, (medio-i), (b.dim-1), funcSort);
		for (i=0; i<(b.dim-1); i++)
			{
				ind = ((b.dim-2) - i);
				HeapSortSparseVector_Swap (b.vpos, b.vpos+(ind+1), b.vval, b.vval+(ind+1));
				HeapSortSparseVector_Adjust(b.vpos, b.vval, 0, ind, funcSort);
			}
	}

void GatherSparseVector (double *v, int n, ptr_SparseVector b)
	{
		int i, dim = 0, *pi1 = b->vpos;
		double *pv1 = b->vval, *pv2 = v;

		for (i=0; i<n; i++) 
			if (*pv2 != 0.0)
				{ *(pi1++) = i; *(pv1++) = *(pv2++); dim++; }
			else pv2++;
		b->dim = dim;
	}

void GatherZeroSparseVector (double *v, int n, ptr_SparseVector b)
	{
		int i, dim = 0, *pi1 = b->vpos;
		double *pv1 = b->vval, *pv2 = v;

		for (i=0; i<n; i++) 
			if (*pv2 != 0.0)
				{ *(pi1++) = i; *(pv1++) = *pv2; *(pv2++) = 0.0; dim++; }
			else pv2++;
		b->dim = dim;
	}

void ScatterSparseVector (SparseVector b, double *v)
	{
		int i, *pi1 = b.vpos;
		double *pv1 = b.vval;

		for (i=0; i<b.dim; i++) v[*(pi1++)] = *(pv1++);
	}

double DotSparseVector (SparseVector b, double *v)
	{
		int i, *pi1 = b.vpos;
		double res = 0.0, *pv1 = b.vval;

		for (i=0; i<b.dim; i++) res += v[*(pi1++)] * *(pv1++);

		return res;
	}

double DotSparseVectors (SparseVector b1, SparseVector b2)
	{
		int i = 0, j = 0, *pi1 = b1.vpos, *pi2 = b2.vpos;
		double res = 0.0, *pv1 = b1.vval, *pv2 = b2.vval;

		while ((i < b1.dim) && (j < b2.dim))
			if (*pi1 < *pi2) { pi1++; pv1++; i++; }
			else if (*pi1 > *pi2) { pi2++; pv2++; j++; }
			else
				{
					res += *(pv1++) * *(pv2++);
					pi1++; i++; pi2++; j++;
				}

		return res;
	}

void AxpySparseVector (SparseVector b, double alfa, double *v)
	{
		int i, *pi1 = b.vpos;
		double *pv1 = b.vval;

		for (i=0; i<b.dim; i++) v[*(pi1++)] += alfa * *(pv1++);
	}

void ScalSparseVector (SparseVector b, double alfa)
	{
		int i = 0;
		double *pd = b.vval;

		for (i=0; i<b.dim; i++) *(pd++) *= alfa;
	}

void VvecSparseVector (SparseVector b, double *v1, double alfa, double *v2)
	{
		int i = 0, *pi = b.vpos;
		double *pd = b.vval;

		for (i=0; i < b.dim; i++)
			{ v2[*pi] += alfa * *(pd++) * v1[*pi]; pi++; }
	}

void VvecSparseVectors (SparseVector b1, SparseVector b2, double alfa, double *v)
	{
		int i = 0, j = 0, *pi1 = b1.vpos, *pi2 = b2.vpos;
		double *pd1 = b1.vval, *pd2 = b2.vval;

		while ((i < b1.dim) && (j < b2.dim))
			if (*pi1 == *pi2)
				{ v[*(pi1++)] += alfa * *(pd1++) * *(pd2++); i++; j++; pi2++; }
			else if (*pi1 < *pi2) { i++; pi1++; pd1++; }
			else { j++; pi2++; pd2++; }
	}

void AddSparseVectors (SparseVector b1, SparseVector b2, ptr_SparseVector b3)
	{
		int i1 = 0, i2 = 0, i3 = 0;
		int n1 = b1.dim, n2 = b2.dim;
		int *pi1 = b1.vpos, *pi2 = b2.vpos, *pi3 = b3->vpos;
		double *pv1 = b1.vval, *pv2 = b2.vval, *pv3 = b3->vval;
	
		if (pi1[n1-1] < *pi2)
			{
				CopySparseVectors (b1, b3); b3->vpos += n1; ; b3->vval += n1; 
				CopySparseVectors (b2, b3); b3->vpos -= n1; ; b3->vval -= n1; 
				b3->dim = (n1+n2);
			}
    else if (pi2[n2-1] < *pi1)
      {
        CopySparseVectors (b2, b3); b3->vpos += n2; ; b3->vval += n2;
        CopySparseVectors (b1, b3); b3->vpos -= n2; ; b3->vval -= n2; 
        b3->dim = (n1+n2);
      }
		else
			{
				while ((i1 < n1) && (i2 < n2))
					if (*pi1 == *pi2)
						{
							*pv3 = *(pv1++) + *(pv2++);
							if (*pv3 != 0.0) { *(pi3++) = *pi1; i3++; pv3++; }
							i1++; pi1++; i2++; pi2++;
						}
					else if (*pi1 < *pi2)
						{ *(pv3++) = *(pv1++); *(pi3++) = *(pi1++); i1++; i3++; }
					else
						{ *(pv3++) = *(pv2++); *(pi3++) = *(pi2++); i2++; i3++; }
				if (i1 < n1)
					{
						CopyInts (pi1, pi3, (n1-i1));
						CopyDoubles (pv1, pv3, (n1-i1));
						i3 += (n1-i1); pi3 += (n1-i1); pv3 += (n1-i1);
					}
				else if (i2 < n2)
					{
						CopyInts (pi2, pi3, (n2-i2));
						CopyDoubles (pv2, pv3, (n2-i2));
						i3 += (n2-i2); pi3 += (n2-i2); pv3 += (n2-i2);
					}
				b3->dim = i3;
			}
	}

void JoinSparseVectors (SparseVector b1, SparseVector b2, ptr_SparseVector b3)
	{
		int i1 = 0, i2 = 0, i3 = 0;
		int n1 = b1.dim, n2 = b2.dim;
		int *pi1 = b1.vpos, *pi2 = b2.vpos, *pi3 = b3->vpos;
		double *pv1 = b1.vval, *pv2 = b2.vval, *pv3 = b3->vval;
	
		if (pi1[n1-1] < *pi2)
			{
				CopySparseVectors (b1, b3); b3->vpos += n1; ; b3->vval += n1; 
				CopySparseVectors (b2, b3); b3->vpos -= n1; ; b3->vval -= n1; 
				b3->dim = (n1+n2);
			}
		else
			{
				while ((i1 < n1) && (i2 < n2))
					if (*pi1 == *pi2)
						{ printf ("ERRROOOOR (JoinSparseVector) \n"); exit (-1); }
					else if (*pi1 < *pi2)
						{ *(pv3++) = *(pv1++); *(pi3++) = *(pi1++); i1++; i3++; }
					else
						{ *(pv3++) = *(pv2++); *(pi3++) = *(pi2++); i2++; i3++; }
				if (i1 < n1)
					{
						CopyInts (pi1, pi3, (n1-i1));
						CopyDoubles (pv1, pv3, (n1-i1));
						i3 += (n1-i1); pi3 += (n1-i1); pv3 += (n1-i1);
					}
				else if (i2 < n2)
					{
            CopyInts (pi2, pi3, (n2-i2));
            CopyDoubles (pv2, pv3, (n2-i2));
						i3 += (n2-i2); pi3 += (n2-i2); pv3 += (n2-i2);
					}
				b3->dim = i3;
			}
	}
