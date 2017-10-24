#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LongVectors.h"

void CreateLongs (long **vlng, int dim)
	{
		if ((*vlng = (long *) malloc (sizeof(long)*dim)) == NULL)
			{ printf ("Memory Error (CreateLongs(%d))\n", dim); exit (1); }
	}

void ReallocLongs (long **vlng, int dim)
	{
		if ((*vlng = (long *) realloc (*vlng, sizeof(long)*dim)) == NULL)
			{ printf ("Memory Error (ReallocLongs(%d))\n", dim); exit (1); }
	}

void InitLongs (long *vlng, int dim, long frst, long incr) 
	{
		int i;
    long *p1 = vlng, num = frst;

		for (i=0; i<dim; i++) 
			{ *(p1++) = num; num += incr; }
	}
		
void CopyLongs (long *src, long *dst, int dim)
	{ 
		memmove (dst, src, sizeof(long) * dim);
	}

void CopyLongsInts (long *src, int *dst, int dim)
	{ 
		int i, *pi = dst;
		long *pl = src;

		for (i=0; i<dim; i++) *(pi++) = (int) *(pl++);
	}

void CopyIntsLongs (int *src, long *dst, int dim)
	{ 
		int i, *pi = src;
		long *pl = dst;

		for (i=0; i<dim; i++) *(pl++) = (long) *(pi++);
	}

void CopyShiftLongs (long *src, long *dst, int dim, long shft)
	{
		int i;
		long *p1 = src, *p2 = dst;

		if (shft == 0)
			CopyLongs (src, dst, dim);
		else
			for (i=0; i<dim; i++)
				*(p2++) = *(p1++) + shft;
	}
	
void GetLongFromString (char *string, long *pnum, int dimC, int shft)
	{
		int j = 0, neg = 0;
		long num = 0;
		char *pchar = string;

		while ((j < dimC) && ((*pchar < '0') || (*pchar > '9')) &&
           (*pchar != '+') && (*pchar != '-')) { j++; pchar++; }
		if (j < dimC)
			{
				if ((*pchar == '+') || (*pchar == '-'))
					{ neg = (*pchar == '-'); j++; pchar++; }
				while ((j < dimC) && (*pchar >= '0') && (*pchar <= '9'))
					{ num = num * 10 + (*pchar - 48); j++; pchar++; }
			}
		if (neg) num = -num;
		*pnum = num + shft; 
	}

void GetLongsFromString (char *string, long *vlng, int dimN, int dimC, int shft)
	{
		int i;
		long *plong = vlng;
		char *pchar = string;

		for (i=0; i<dimN; i++)
			{ GetLongFromString (pchar, (plong++), dimC, shft); pchar += dimC; }
	}

void TransformLengthtoHeaderL (long *vlng, int dim)
	{
		int i;
		long *pl = vlng; 
    for (i=0; i<dim; i++) { *(pl+1) += *pl; pl++; }
	}

void TransformHeadertoLengthL (long *vlng, int dim)
	{
		int i;
		long *pl = vlng+dim; 
    for (i=dim; i>0; i--) { *(pl) -= *(pl-1); pl--; }
	}

void ComputeHeaderfromLengthL (long *len, long *head, int dim)
	{
		int i;
		long *pl1 = len, *pl2 = head; 
    for (i=0; i<dim; i++) { *(pl2+1) = (*pl2) +(*(pl1++)); pl2++; }
	}

void ComputeLengthfromHeaderL (long *head, long *len, int dim)
	{
		int i;
		long *pl1 = head, *pl2 = len; 
    for (i=0; i<dim; i++) { *(pl2++) = (*(pl1+1)) -(*pl1); pl1++; }
	}

int HeapSortPermLongs_IncrPos (long n1, long n2, int p1, int p2)
	{ return ((n1 < n2) || ((n1 == n2) && (p1 < p2))); }

int HeapSortPermLongs_DecrPos (long n1, long n2, int p1, int p2)
	{ return ((n1 > n2) || ((n1 == n2) && (p1 > p2))); }

int HeapSortPermLongs_IncrVal (long n1, long n2, int p1, int p2)
	{ return (p1 < p2); }

int HeapSortPermLongs_DecrVal (long n1, long n2, int p1, int p2)
	{ return (p1 > p2); }

void HeapSortPermLongs_Swap (long *n1, long *n2, int *p1, int *p2)
	{
		long n; 
		int p;
		n = *n1; *n1 = *n2; *n2 = n;
		p = *p1; *p1 = *p2; *p2 = p;
	}

void HeapSortPermLongs_Adjust (long *vlng, int *prm, int inic, int final, Longs_SortPerm_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vlng[j], vlng[j+1], prm[j], prm[j+1]));
				_bool = funcSort (vlng[k], vlng[j], prm[k], prm[j]);
				if (_bool)
					{ 
						HeapSortPermLongs_Swap (vlng+k, vlng+j, prm+k, prm+j); 
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortPermLongs (long *vlng, int *prm, int dim, Longs_SortPerm_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortPermLongs_Adjust (vlng, prm, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortPermLongs_Swap (vlng, vlng+(ind+1), prm, prm+(ind+1));
				HeapSortPermLongs_Adjust(vlng, prm, 0, ind, funcSort);
			}
	}

int HeapSortLongs_IncrPos (long n1, long n2)
	{ return (n1 < n2); }

int HeapSortLongs_DecrPos (long n1, long n2)
	{ return (n1 > n2); }

void HeapSortLongs_Swap (long *n1, long *n2)
	{
		long n;
		n = *n1; *n1 = *n2; *n2 = n;
	}

void HeapSortLongs_Adjust (long *vlng, int inic, int final, Longs_Sort_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vlng[j], vlng[j+1]));
				_bool = funcSort (vlng[k], vlng[j]);
				if (_bool)
					{ 
						HeapSortLongs_Swap (vlng+k, vlng+j);
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortLongs (long *vlng, int dim, Longs_Sort_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortLongs_Adjust (vlng, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortLongs_Swap (vlng, vlng+(ind+1));
				HeapSortLongs_Adjust(vlng, 0, ind, funcSort);
			}
	}

int FindSortLong (long num, long *vlng, int dim, Longs_Sort_func funcSort)
	{
		int _bool = 1, izq = 0, der = dim - 1, cnt = 0; 

		while ((_bool) && (izq <= der))
			{
				cnt = (izq + der) / 2;
				if (num == vlng[cnt])
					_bool = 0;
				else if (funcSort (num, vlng[cnt])) der = cnt - 1;
				else izq = cnt + 1;
			}
		if (_bool) cnt = -1;
		else if (cnt > 0)
			{
				if (funcSort (0,1))
					{ while ((cnt >= 0) && (vlng[cnt] == num)) cnt--; cnt++; }
				else
					{ while ((cnt < dim) && (vlng[cnt] == num)) cnt++; cnt--; }
			}
		return cnt;
	}
	
long AddLongs (long *vlng, int dim) 
	{
		int i;
		long *pl = vlng, aux = 0;
		for (i=0; i<dim; i++) { aux += *pl; pl++; }
		return aux;
	}

long AddPermuteLongs (long *vlng, int *perm, int dim) 
	{
		int i, *pi = perm;
		long aux = 0;
		
		for (i=0; i<dim; i++) { aux += vlng[*pi]; pi++; }
		return aux;
	}

long MaxLongs (long *vlng, int dim) 
	{
		int i;
		long *pl = vlng, aux = *pl;

		for (i=1; i<dim; i++) { if (*pl > aux) aux = *pl; pl++; }
		return aux;
	}

void PermuteLongs (long *vlng, int *perm, int dim)
	{
		int i;
		long *pl = vlng;
		
		for (i=0; i<dim; i++) { *pl = perm[*pl]; pl++; }
	}

void CopyPermuteLongs (long *src, int *perm, long *dst, int dim)
	{
		int i, *pp = perm;
		long *pi1 = src, *pi2 = dst;
		
		for (i=0; i<dim; i++) { *pi2 = pi1[*pp]; pp++; pi2++; }
	}

void CopyInvPermuteLongs (long *src, long *dst, int *perm, int dim)
	{
		int i, *pp = perm;
		long *pl1 = src, *pl2 = dst;
		
		for (i=0; i<dim; i++) { pl2[*pp] = *pl1; pp++; pl1++; }
	}

void PrintLongs (long *vlng, int dim)
	{
		int i;
		long *pl = vlng;

		for (i=0; i<dim; i++) printf ("%ld ", *(pl++));
		printf ("\n");
	}

void FPrintLongs (FILE * f, long *vlng, int dim)
	{
		int i;
		long *pl = vlng;

		for (i=0; i<dim; i++) fprintf (f, "%ld ", *(pl++));
		fprintf (f, "\n");
	}

void PrintFLongs (long *vlng, int dim, int f1, int f2) {
	char formato[10];
	int i;
	long *pf = vlng;

	if (f2 > 0)
		sprintf (formato, "%%%d.%dld", f1,f2);
	else
		sprintf (formato, "%%%dld", f1);

	for (i=0; i<dim; i++) printf (formato, *(pf++));
	printf ("\n");
}

int ReadLongs (char *file, long **vlng) {
	int n = -1;
	long num;
	FILE *f = NULL;

	f = fopen (file, "r");
	if (f != NULL) {
		n = 0;
		while (fscanf (f, "%ld", &num) != EOF) n++;
		fclose (f);

		CreateLongs (vlng, n); n = 0;
		f = fopen (file, "r"); 
		while (fscanf (f, "%ld", &num) != EOF) 
			(*vlng)[n++] = num;
		fclose (f);
	}

	return n;
}

void WriteLongs (char *filename, long *vlng, int dim) {
	FILE *f = NULL;

	printf ("Print %s of size %d\n", filename, dim);
	f = fopen (filename, "w");
	FPrintLongs (f, vlng, dim);
	fclose (f);
}

void RemoveLongs (long **vlng)
	{ if (*vlng != NULL) free (*vlng); *vlng = NULL; }

/***********************************************************/

void CreateVectorPtrLongs (matLongs *mlng, int numR)
	{
		if ((*mlng = (matLongs) malloc (sizeof(long *) * numR)) == NULL)
			{ printf ("Memory Error (CreateMatrixLongs(%d))\n", numR); exit (1); }
	}

void CopyVectorPtrLongs (matLongs msrc, matLongs mdst, int dim)
	{ memmove (mdst, msrc, sizeof(long *) * dim); }

void DesplVectorPtrLongs (matLongs mlng, int dim, int dspl)
	{ 
		int i; 
		for (i=0; i<dim; i++) mlng[i] += dspl; 
	}

void RemoveVectorPtrLongs (matLongs *mlng)
	{
		if (*mlng != NULL) free (*mlng); 
		*mlng = NULL;
	}

void CreateMatrixLongs (matLongs *mlng, int numR, int numC)
	{
		int i;

		CreateVectorPtrLongs (mlng, numR);
		CreateLongs (*mlng, numR*numC);
		for (i=1; i<numR; i++) (*mlng)[i] = (*mlng)[i-1] + numC;
	}

void RemoveMatrixLongs (matLongs *mlng)
	{
		if (*mlng != NULL) RemoveLongs (*mlng); 
		RemoveVectorPtrLongs (mlng);
	}

