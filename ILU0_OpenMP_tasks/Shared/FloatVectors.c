#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FloatVectors.h"

void CreateFloats (float **vflt, int dim)
	{
		if ((*vflt = (float *) malloc (sizeof(float)*dim)) == NULL)
			{ printf ("Memory Error (CreateFloats(%d))\n", dim); exit (1); }
	}

void ReallocFloats (float **vftl, int dim)
	{
		if ((*vftl = (float *) realloc (*vftl, sizeof(float)*dim)) == NULL)
			{ printf ("Memory Error (ReallocFloats(%d))\n", dim); exit (1); }
	}

void InitFloats (float *vflt, int dim, float frst, float incr) 
	{
		int i; 
		float *pf = vflt, num = frst;

		for (i=0; i<dim; i++) 
			{ *(pf++) = num; num += incr; }
	}
		
void CopyFloats (float *src, float *dst, int dim)
	{ 
		memmove (dst, src, sizeof(float) * dim);
	}

void ScaleFloats (float *vflt, float scal, int dim) 
	{
		int i; 
		float *pf = vflt;

		for (i=0; i<dim; i++) 
			*(pf++) *= scal;
	}

void AxpyFloats (float alfa, float *vflt1, float *vflt2, int dim) 
	{
		int i; 
		float *pf1 = vflt1, *pf2 = vflt2;

		for (i=0; i<dim; i++) 
			*(pf2++) += (*(pf1++) * alfa);
	}

void XpayFloats (float *vflt1, float alfa, float *vflt2, int dim) 
	{
		int i; 
		float *pf1 = vflt1, *pf2 = vflt2;

		for (i=0; i<dim; i++) 
			{ *pf2 = ((*pf2) * alfa) + (*(pf1++)); pf2++; }
	}

float DotFloats (float *vflt1, float *vflt2, int dim) 
	{
		int i; 
		float *pf1 = vflt1, *pf2 = vflt2, res = 0.0;

		for (i=0; i<dim; i++) 
			res += (*(pf1++)) * (*(pf2++));

		return res;
	}

void GetFloatFromString (char *string, float *pflt, int dimC)
	{
		int j, k, exp, neg;
		float num, frac;
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
		*pflt = num; 
	}

void GetFloatsFromString (char *string, float *vflt, int dimN, int dimC)
	{
		int i;
		float *paux = vflt;
		char *pchar = string;

		for (i=0; i<dimN; i++)
			{ GetFloatFromString (pchar, (paux++), dimC); pchar += dimC; }
	}

int HeapSortPermFloats_IncrPos (float n1, float n2, int p1, int p2)
	{ return ((n1 < n2) || ((n1 == n2) && (p1 < p2))); }

int HeapSortPermFloats_DecrPos (float n1, float n2, int p1, int p2)
	{ return ((n1 > n2) || ((n1 == n2) && (p1 > p2))); }

int HeapSortPermFloats_IncrVal (float n1, float n2, int p1, int p2)
	{ return (p1 < p2); }

int HeapSortPermFloats_DecrVal (float n1, float n2, int p1, int p2)
	{ return (p1 > p2); }

void HeapSortPermFloats_Swap (float *n1, float *n2, int *p1, int *p2)
	{
		float n;
		int p;
		n = *n1; *n1 = *n2; *n2 = n;
		p = *p1; *p1 = *p2; *p2 = p;
	}

void HeapSortPermFloats_Adjust (float *vflt, int *prm, int inic, int final, Floats_SortPerm_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vflt[j], vflt[j+1], prm[j], prm[j+1]));
				_bool = funcSort (vflt[k], vflt[j], prm[k], prm[j]);
				if (_bool)
					{ 
						HeapSortPermFloats_Swap (vflt+k, vflt+j, prm+k, prm+j); 
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortPermFloats (float *vflt, int *prm, int dim, Floats_SortPerm_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortPermFloats_Adjust (vflt, prm, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortPermFloats_Swap (vflt, vflt+(ind+1), prm, prm+(ind+1));
				HeapSortPermFloats_Adjust(vflt, prm, 0, ind, funcSort);
			}
	}

int HeapSortFloats_IncrPos (float n1, float n2)
	{ return (n1 < n2); }

int HeapSortFloats_DecrPos (float n1, float n2)
	{ return (n1 > n2); }

void HeapSortFloats_Swap (float *n1, float *n2)
	{
		float n;
		n = *n1; *n1 = *n2; *n2 = n;
	}

void HeapSortFloats_Adjust (float *vflt, int inic, int final, Floats_Sort_func funcSort)
	{
		int j, k, _bool;
		
		k = inic; j = (k << 1) + 1; _bool = 1;
		while ((j <= final) && _bool)
			{
				j += ((j < final) && funcSort (vflt[j], vflt[j+1]));
				_bool = funcSort (vflt[k], vflt[j]);
				if (_bool)
					{ 
						HeapSortFloats_Swap (vflt+k, vflt+j);
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortFloats (float *vflt, int dim, Floats_Sort_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortFloats_Adjust (vflt, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortFloats_Swap (vflt, vflt+(ind+1));
				HeapSortFloats_Adjust(vflt, 0, ind, funcSort);
			}
	}

int FindSortFloat (float num, float *vflt, int dim, Floats_Sort_func funcSort)
	{
		int _bool = 1, izq = 0, der = dim - 1, cnt = 0; 

		while ((_bool) && (izq <= der))
			{
				cnt = (izq + der) / 2;
				if (num == vflt[cnt])
					_bool = 0;
				else if (funcSort (num, vflt[cnt])) der = cnt - 1;
				else izq = cnt + 1;
			}
		if (_bool) cnt = -1;
		else if (cnt > 0)
			{
				if (funcSort (0,1))
					{ while ((cnt >= 0) && (vflt[cnt] == num)) cnt--; cnt++; }
				else
					{ while ((cnt < dim) && (vflt[cnt] == num)) cnt++; cnt--; }
			}
		return cnt;
	}
	
void VdivFloats (float alfa, float *src, float *dst, int dim)
  {
    int i;

    for (i=0; i<dim; i++)
      { *dst = (alfa * *dst / *(src++)); dst++; } 
  }

void VvecFloats (float alfa, float *src1, float *src2, float beta, float *dst, int dim)
  {
    int i;

    for (i=0; i<dim; i++)
      { *dst = (beta * *dst) + (alfa * *(src1++) * *(src2++)); dst++; }
  }

float AddFloats (float *vflt, int dim) 
	{
		int i;
		float *pf = vflt, aux = 0;
		for (i=0; i<dim; i++) { aux += *pf; pf++; }
		return aux;
	}

float AddPermuteFloats (float *vflt, int *perm, int dim) 
	{
		int i, *pi = perm;
		float aux = 0;
		
		for (i=0; i<dim; i++) { aux += vflt[*pi]; pi++; }
		return aux;
	}

float MaxFloats (float *vflt, int dim) 
	{
		int i;
		float *pf = vflt, aux = *pf;

		for (i=1; i<dim; i++) { if (*pf > aux) aux = *pf; pf++; }
		return aux;
	}

void CopyPermuteFloats (float *src, int *perm, float *dst, int dim)
	{
		int i, *pp = perm;
		float *pf1 = src, *pf2 = dst;
		
		for (i=0; i<dim; i++) { *pf2 = pf1[*pp]; pp++; pf2++; }
	}

void CopyInvPermuteFloats (float *src, float *dst, int *perm, int dim)
	{
		int i, *pp = perm;
		float *pf1 = src, *pf2 = dst;
		
		for (i=0; i<dim; i++) { pf2[*pp] = *pf1; pp++; pf1++; }
	}

void PrintFloats (float *vflt, int dim)
	{
		int i;
		float *pf = vflt;

		for (i=0; i<dim; i++) printf ("%e ", *(pf++));
		printf ("\n");
	}

void FPrintFloats (FILE *f, float *vflt, int dim)
	{
		int i;
		float *pf = vflt;

		for (i=0; i<dim; i++) fprintf (f, "%20.14g ", *(pf++));
		fprintf (f, "\n");
	}

void PrintFFloats (float *vflt, int dim, int f1, int f2) {
	char formato[10];
	int i;
	float *pf = vflt;

	if (f2 > 0)
		sprintf (formato, "%%%d.%de", f1,f2);
	else
		sprintf (formato, "%%%de", f1);

	for (i=0; i<dim; i++) printf (formato, *(pf++));
	printf ("\n");
}

int ReadFloats (char *file, float **vflt) {
	int n = -1;
	float num;
	FILE *f = NULL;

	f = fopen (file, "r");
	if (f != NULL) {
		n = 0;
		while (fscanf (f, "%f", &num) != EOF) n++;
		fclose (f);
	
		CreateFloats (vflt, n); n = 0;
		f = fopen (file, "r"); 
		while (fscanf (f, "%f", &num) != EOF) 
			(*vflt)[n++] = num;
		fclose (f);
	}

	return n;
}

void WriteFloats (char *filename, float *vflt, int dim) {
	FILE *f = NULL;

	printf ("Print %s of size %d\n", filename, dim);
	f = fopen (filename, "w");
	FPrintFloats (f, vflt, dim);
	fclose (f);
}

void RemoveFloats (float **vflt)
	{ if (*vflt != NULL) free (*vflt); *vflt = NULL; }

/****************************************************************************/

void CreateVectorPtrFloats (matFloats *mflt, int numR)
	{
		if ((*mflt = (matFloats) malloc (sizeof(float *) * numR)) == NULL)
			{ printf ("Memory Error (CreateMatrixFloats(%d))\n", numR); exit (1); }
	}

void CopyVectorPtrFloats (matFloats msrc, matFloats mdst, int dim)
	{ memmove (mdst, msrc, sizeof(float *) * dim); }

void DesplVectorPtrFloats (matFloats mflt, int dim, int dspl)
	{ 
		int i; 
		for (i=0; i<dim; i++) mflt[i] += dspl; 
	}

void RemoveVectorPtrFloats (matFloats *mflt)
	{
		if (*mflt != NULL) free (*mflt); 
		*mflt = NULL;
	}

void CreateMatrixFloats (matFloats *mflt, int numR, int numC)
	{
		int i;

		CreateVectorPtrFloats (mflt, numR);
		CreateFloats (*mflt, numR*numC);
		for (i=1; i<numR; i++) (*mflt)[i] = (*mflt)[i-1] + numC;
	}

void RemoveMatrixFloats (matFloats *mflt)
	{
		if (*mflt != NULL) RemoveFloats (*mflt); 
		RemoveVectorPtrFloats (mflt);
	}

