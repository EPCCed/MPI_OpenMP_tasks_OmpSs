#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "IntVectors.h"

void CreateInts (int **vint, int dim)
	{
		if ((*vint = (int *) malloc (sizeof(int)*dim)) == NULL)
			{ printf ("Memory Error (CreateInts(%d))\n", dim); exit (1); }
	}

void ReallocInts (int **vint, int dim)
	{
		if ((*vint = (int *) realloc (*vint, sizeof(int)*dim)) == NULL)
			{ printf ("Memory Error (ReallocInts(%d))\n", dim); exit (1); }
	}

void InitInts (int *vint, int dim, int frst, int incr) 
	{
		int i, *p1 = vint, num = frst;

		for (i=0; i<dim; i++) 
			{ *(p1++) = num; num += incr; }
	}
		
void CopyInts (int *src, int *dst, int dim)
	{ 
		memmove (dst, src, sizeof(int) * dim);
	}

void CopyShiftInts (int *src, int *dst, int dim, int shft)
	{
		int i, *p1 = src, *p2 = dst;

		if (shft == 0)
			CopyInts (src, dst, dim);
		else
			for (i=0; i<dim; i++)
				*(p2++) = *(p1++) + shft;
	}
	
void CopyShiftRangsInts (int *src, int *dst, int dim1,int *rangs,  int *shfts, int dim2) 
	{
		int i, *p1 = src, *p2 = dst, j, *q1 = rangs, *q2 = shfts, *q3 = NULL;

		while ((dim2 > 0) && (*q2 == 0)) { dim2--; q1++; q2++; }

		if (dim2 == 0)
			CopyInts (src, dst, dim1);
		else 
			{
				for (i=0; i<dim1; i++)
					{
						j = 0; q3 = q1;
						while ((j < dim2) && (*q3 <= *p1)) { j++; q3++; }
						if (j > 0)
							*(p2++) = *(p1++) + q2[--j];
						else
							*(p2++) = *(p1++);
					}
			}
	}
	
void GetIntFromString (char *string, int *pnum, int dimC, int shft)
	{
		int j = 0, num = 0, neg = 0;
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

void GetIntsFromString (char *string, int *vint, int dimN, int dimC, int shft)
	{
		int i, *pint = vint;
		char *pchar = string;

		for (i=0; i<dimN; i++)
			{ GetIntFromString (pchar, (pint++), dimC, shft); pchar += dimC; }
	}

void GetFormatsFromString (char *string, int *vint, int dimN, int dimC)
	{
		int i, k = 0;
		int *pint = vint;
		char *pchar = string, *pch = NULL, c = ' ', c2;

		for (i=0; i<dimN; i++)
			{ 
				pch = pchar; 
				while (*pch == ' ') pch++;
				sscanf (pch, "(%i%c", pint, &c);
				if ((c == 'P') || (c == 'p')) 
					{ 
						sscanf (pch, "(%i%c%i%c%i.%i)", &k, &c2, pint, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else if ((c == 'E') || (c == 'e') || (c == 'D') || (c == 'd')	||
								 (c == 'F') || (c == 'f') || (c == 'G') || (c == 'g'))
					{ 
						sscanf (pch, "(%i%c%i.%i)", pint, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else 
					{ sscanf (pch, "(%i%c%i)", pint, &c, pint+1); pint += 2; }
				pchar += dimC;
			}
	}

void GetFormatsFromString2 (char *string, int *vint, int dimN, int dimC)
	{
		int i, k = 0;
		int *pint = vint;
		char *pchar = string, *pch = NULL, c = ' ', c2;

		for (i=0; i<dimN; i++)
			{ 
				pch = pchar; 
				while (*pch == ' ') pch++;
				sscanf (pch, "(%i%c", pint, &c);
				if (c == '(')
					{ 
						sscanf (pch, "(%i(%i%c,%c%i.%i))", pint, &k, &c2, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else if ((c == 'P') || (c == 'p')) 
					{ 
						sscanf (pch, "(%i%c%i%c%i.%i)", &k, &c2, pint, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else if ((c == 'E') || (c == 'e') || (c == 'D') || (c == 'd')	||
								 (c == 'F') || (c == 'f') || (c == 'G') || (c == 'g'))
					{ 
						sscanf (pch, "(%i%c%i.%i)", pint, &c, pint+1, pint+2); 
						pint += 3; 
					}
				else if ((c == 'X') || (c == 'x')) 
					{ sscanf (pch, "(%i%c,%i%c%i)", pint, &c2, pint+1, &c, pint+2); pint += 3; }
				else 
					{ *pint = 0; sscanf (pch, "(%i%c%i)", pint+1, &c, pint+2); pint += 3; }
				pchar += dimC;
			}
	}

void TransformLengthtoHeader (int *vint, int dim)
	{
		int i, *pi = vint; 
    for (i=0; i<dim; i++) { *(pi+1) += *pi; pi++; }
	}

void TransformHeadertoLength (int *vint, int dim)
	{
		int i, *pi = vint+dim; 
    for (i=dim; i>0; i--) { *(pi) -= *(pi-1); pi--; }
	}

void ComputeHeaderfromLength (int *len, int *head, int dim)
	{
		int i, *pi1 = len, *pi2 = head; 
    for (i=0; i<dim; i++) { *(pi2+1) = (*pi2) +(*(pi1++)); pi2++; }
	}

void ComputeLengthfromHeader (int *head, int *len, int dim)
	{
		int i, *pi1 = head, *pi2 = len; 
    for (i=0; i<dim; i++) { *(pi2++) = (*(pi1+1)) -(*pi1); pi1++; }
	}

int HeapSortPermInts_IncrPos (int n1, int n2, int p1, int p2)
	{ return ((n1 < n2) || ((n1 == n2) && (p1 < p2))); }

int HeapSortPermInts_DecrPos (int n1, int n2, int p1, int p2)
	{ return ((n1 > n2) || ((n1 == n2) && (p1 > p2))); }

int HeapSortPermInts_IncrVal (int n1, int n2, int p1, int p2)
	{ return (p1 < p2); }

int HeapSortPermInts_DecrVal (int n1, int n2, int p1, int p2)
	{ return (p1 > p2); }

void HeapSortPermInts_Swap (int *n1, int *n2, int *p1, int *p2)
	{
		int n, p;
		n = *n1; *n1 = *n2; *n2 = n;
		p = *p1; *p1 = *p2; *p2 = p;
	}

void HeapSortPermInts_Adjust (int *vint, int *prm, int inic, int final, Ints_SortPerm_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vint[j], vint[j+1], prm[j], prm[j+1]));
				_bool = funcSort (vint[k], vint[j], prm[k], prm[j]);
				if (_bool)
					{ 
						HeapSortPermInts_Swap (vint+k, vint+j, prm+k, prm+j); 
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortPermInts (int *vint, int *prm, int dim, Ints_SortPerm_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortPermInts_Adjust (vint, prm, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortPermInts_Swap (vint, vint+(ind+1), prm, prm+(ind+1));
				HeapSortPermInts_Adjust(vint, prm, 0, ind, funcSort);
			}
	}

int HeapSortInts_IncrPos (int n1, int n2)
	{ return (n1 < n2); }

int HeapSortInts_DecrPos (int n1, int n2)
	{ return (n1 > n2); }

void HeapSortInts_Swap (int *n1, int *n2)
	{
		int n;
		n = *n1; *n1 = *n2; *n2 = n;
	}

void HeapSortInts_Adjust (int *vint, int inic, int final, Ints_Sort_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vint[j], vint[j+1]));
				_bool = funcSort (vint[k], vint[j]);
				if (_bool)
					{ 
						HeapSortInts_Swap (vint+k, vint+j);
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortInts (int *vint, int dim, Ints_Sort_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortInts_Adjust (vint, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortInts_Swap (vint, vint+(ind+1));
				HeapSortInts_Adjust(vint, 0, ind, funcSort);
			}
	}

int FindSortInt (int num, int *vint, int dim, Ints_Sort_func funcSort)
	{
		int _bool = 1, izq = 0, der = dim - 1, cnt = 0; 

		while ((_bool) && (izq <= der))
			{
				cnt = (izq + der) / 2;
				if (num == vint[cnt])
					_bool = 0;
				else if (funcSort (num, vint[cnt])) der = cnt - 1;
				else izq = cnt + 1;
			}
		if (_bool) cnt = -1;
		else if (cnt > 0)
			{
				if (funcSort (0,1))
					{ while ((cnt >= 0) && (vint[cnt] == num)) cnt--; cnt++; }
				else
					{ while ((cnt < dim) && (vint[cnt] == num)) cnt++; cnt--; }
			}
		return cnt;
	}
	
int AddInts (int *vint, int dim) 
	{
		int i, *pi = vint, aux = 0;
		for (i=0; i<dim; i++) { aux += *pi; pi++; }
		return aux;
	}

int AddPermuteInts (int *vint, int *perm, int index, int dim) 
	{
		int i, *pi1 = perm, *pi2 = vint - index, aux = 0;
		
		for (i=0; i<dim; i++) { aux += pi2[*pi1]; pi1++; }
		return aux;
	}

int MaxInts (int *vint, int dim) 
	{
		int i, *pi = vint, aux = *pi;

		for (i=1; i<dim; i++) { if (*pi > aux) aux = *pi; pi++; }
		return aux;
	}

void PermuteInts (int *vint, int *perm, int index, int dim)
	{
		int i, *pi1 = vint, *pi2 = perm - index;
		
		for (i=0; i<dim; i++) { *pi1 = pi2[*pi1]; pi1++; }
	}

void CopyPermuteInts (int *src, int *perm, int index, int *dst, int dim)
	{
		int i, *pi1 = src - index, *pp = perm, *pi2 = dst;
		
		for (i=0; i<dim; i++) { *pi2 = pi1[*pp]; pp++; pi2++; }
	}

void CopyInvPermuteInts (int *src, int *dst, int *perm, int index, int dim)
	{
		int i, *pp = perm, *pi1 = src, *pi2 = dst - index;
		
		for (i=0; i<dim; i++) { pi2[*pp] = *pi1; pp++; pi1++; }
	}

void ComputeInvPermutation (int *perm, int *iperm, int index, int dim)
	{
		int i, *pi1 = perm, *pi2 = iperm - index;
		
		for (i=index; i<(dim+index); i++) { pi2[*(pi1++)] = i; }
	}

void PrintInts (int *vint, int dim)
	{
		int i, *pi = vint;

		for (i=0; i<dim; i++) printf ("%d ", *(pi++));
		printf ("\n");
	}

void PrintSeqInts (int inic, int dim)
	{
		int i, num = inic;

		for (i=0; i<dim; i++) printf ("%d ", num++);
		printf ("\n");
	}

void FPrintInts (FILE * f, int *vint, int dim)
	{
		int i, *pi = vint;

		for (i=0; i<dim; i++) fprintf (f, "%d ", *(pi++));
		fprintf (f, "\n");
	}

void FPrintSeqInts (FILE * f, int inic, int dim)
	{
		int i, num = inic;

		for (i=0; i<dim; i++) fprintf (f, "%d ", num++);
		fprintf (f, "\n");
	}

void PrintFInts (int *vint, int dim, int f1, int f2) {
	char formato[10];
	int i, *pi = vint;

	if (f1 <= 0)
		sprintf (formato, "%%d ");
	else if (f2 <= 0)
		sprintf (formato, "%%%dd", f1);
	else
		sprintf (formato, "%%%d.%dd", f1,f2);

	for (i=0; i<dim; i++) printf (formato, *(pi++));
	printf ("\n");
}

void PrintFSeqInts (int inic, int dim, int f1, int f2) {
	char formato[10];
	int i, num = inic;

	if (f1 <= 0)
		sprintf (formato, "%%d ");
	else if (f2 <= 0)
		sprintf (formato, "%%%dd", f1);
	else
		sprintf (formato, "%%%d.%dd", f1,f2);

	for (i=0; i<dim; i++) printf (formato, num++);
	printf ("\n");
}

void FPrintFInts (FILE *f, int *vint, int dim, int f1, int f2) {
	char formato[10];
	int i, j, num, *pi = vint;

	if (f1 <= 0)
		sprintf (formato, "%%d ");
	else if (f2 <= 0)
		sprintf (formato, "%%%dd", f1);
	else
		sprintf (formato, "%%%d.%dd", f1,f2);

	num = (f1 <=0)? 15: (80 / f1); j = num; 
	for (i=0; i<dim; i++) {
		fprintf (f, formato, *(pi++));
		if ((--j) == 0) {
			fprintf (f, "\n"); j = num;
		}
	}		
	if (j < num) fprintf (f, "\n");
}

int ReadInts (char *file, int **vint) {
	int n = -1, num;
	FILE *f = NULL;

	f = fopen (file, "r");
	if (f != NULL) {
		n = 0;
		while (fscanf (f, "%d", &num) != EOF) n++;
		fclose (f);

		CreateInts (vint, n); n = 0;
		f = fopen (file, "r"); 
		while (fscanf (f, "%d", &num) != EOF) 
			(*vint)[n++] = num;
		fclose (f);
	}

	return n;
}

void WriteInts (char *filename, int *vint, int dim) {
	FILE *f = NULL;

	printf ("Print %s of size %d\n", filename, dim);
	f = fopen (filename, "w");
	FPrintInts (f, vint, dim);
	fclose (f);
}

void WriteFInts (char *filename, int *vint, int dim, int f1, int f2) {
	FILE *f = NULL;

	printf ("Print %s of size %d\n", filename, dim);
	f = fopen (filename, "w");
	FPrintFInts (f, vint, dim, f1, f2);
	fclose (f);
}

void RemoveInts (int **vint)
	{ if (*vint != NULL) free (*vint); *vint = NULL; }

/***********************************************************/

void CreateVectorPtrInts (matInts *mint, int numR)
	{
		if ((*mint = (matInts) malloc (sizeof(int *) * numR)) == NULL)
			{ printf ("Memory Error (CreateMatrixInts(%d))\n", numR); exit (1); }
	}

void CopyVectorPtrInts (matInts msrc, matInts mdst, int dim)
	{ memmove (mdst, msrc, sizeof(int *) * dim); }

void DesplVectorPtrInts (matInts mint, int dim, int dspl)
	{ 
		int i; 
		for (i=0; i<dim; i++) mint[i] += dspl; 
	}

void RemoveVectorPtrInts (matInts *mint)
	{
		if (*mint != NULL) free (*mint); 
		*mint = NULL;
	}

void CreateMatrixInts (matInts *mint, int numR, int numC)
	{
		int i;

		CreateVectorPtrInts (mint, numR);
		CreateInts (*mint, numR*numC);
		for (i=1; i<numR; i++) (*mint)[i] = (*mint)[i-1] + numC;
	}

void RemoveMatrixInts (matInts *mint)
	{
		if (*mint != NULL) RemoveInts (*mint); 
		RemoveVectorPtrInts (mint);
	}

