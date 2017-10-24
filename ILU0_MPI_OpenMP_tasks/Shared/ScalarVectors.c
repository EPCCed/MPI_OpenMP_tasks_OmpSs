#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ScalarVectors.h"

int CreateInts (int **vint, int num)
	{
		if ((*vint = (int *) malloc (sizeof(int)*num)) == NULL)
			{ printf ("Memory Error (CreateInts(%d))\n", num); exit (1); }
	}

int InitInts (int *vint, int n, int frst, int incr) 
	{
		int i, *p1 = vint, num = frst;

		for (i=0; i<n; i++) 
			{ *(p1++) = num; num += incr; }
	}
		
int CopyInts (int *src, int *dst, int n)
	{ 
		memmove (dst, src, sizeof(int) * n);
	}

int CopyShiftInts (int *src, int *dst, int n, int shft)
	{
		int i, *p1 = src, *p2 = dst;

		if (shft == 0)
			CopyInts (src, dst, n);
		else
			for (i=0; i<n; i++)
				*(p2++) = *(p1++) + shft;
	}
	
void GetIntFromString (char *string, int *pnum, int numC, int shft)
	{
		int j = 0, num = 0, neg = 0;
		char *pchar = string;

		while ((j < numC) && ((*pchar < '0') || (*pchar > '9')) &&
           (*pchar != '+') && (*pchar != '-')) { j++; pchar++; }
		if (j < numC)
			{
				if ((*pchar == '+') || (*pchar == '-'))
					{ neg = (*pchar == '-'); j++; pchar++; }
				while ((j < numC) && (*pchar >= '0') && (*pchar <= '9'))
					{ num = num * 10 + (*pchar - 48); j++; pchar++; }
			}
		if (neg) num = -num;
		*pnum = num + shft; 
	}

void GetIntsFromString (char *string, int *vec, int numN, int numC, int shft)
	{
		int i, *pint = vec;
		char *pchar = string;

		for (i=0; i<numN; i++)
			{ GetIntFromString (pchar, (pint++), numC, shft); pchar += numC; }
	}

void GetFormatsFromString (char *string, int *vec, int numN, int numC)
	{
		int i, k = 0;
		int *pint = vec;
		char *pchar = string, *pch = NULL, c = ' ', c2;

		for (i=0; i<numN; i++)
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
				pchar += numC;
			}
	}

void GetFormatsFromString2 (char *string, int *vec, int numN, int numC)
	{
		int i, k = 0;
		int *pint = vec;
		char *pchar = string, *pch = NULL, c = ' ', c2;

		for (i=0; i<numN; i++)
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
				pchar += numC;
			}
	}

int TransformLengthtoHeader (int *vec, int n)
	{
		int i, *pi = vec; 
    for (i=0; i<n; i++) { *(pi+1) += *pi; pi++; }
	}

int TransformHeadertoLength (int *vec, int n)
	{
		int i, *pi = vec+n; 
    for (i=n; i>0; i--) { *(pi) -= *(pi-1); pi--; }
	}

int ComputeHeaderfromLength (int *len, int *head, int n)
	{
		int i, *pi1 = len, *pi2 = head; 
    for (i=0; i<n; i++) { *(pi2+1) = (*pi2) +(*(pi1++)); pi2++; }
	}

int ComputeLengthfromHeader (int *head, int *len, int n)
	{
		int i, *pi1 = head, *pi2 = len; 
    for (i=0; i<n; i++) { *(pi2++) = (*(pi1+1)) -(*pi1); pi1++; }
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

void HeapSortPermInts_Adjust (int *vec, int *prm, int inic, int final, Ints_SortPerm_func funcSort)
	{
		int j, k, bool;
		
		k = inic; j = (k << 1) + 1; bool = 1;
		while ((j <= final) && bool)
			{
				j += ((j < final) && funcSort (vec[j], vec[j+1], prm[j], prm[j+1]));
				bool = funcSort (vec[k], vec[j], prm[k], prm[j]);
				if (bool)
					{ 
						HeapSortPermInts_Swap (vec+k, vec+j, prm+k, prm+j); 
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortPermInts (int *vec, int *prm, int dim, Ints_SortPerm_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortPermInts_Adjust (vec, prm, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortPermInts_Swap (vec, vec+(ind+1), prm, prm+(ind+1));
				HeapSortPermInts_Adjust(vec, prm, 0, ind, funcSort);
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

void HeapSortInts_Adjust (int *vec, int inic, int final, Ints_Sort_func funcSort)
	{
		int j, k, bool;
		
		k = inic; j = (k << 1) + 1; bool = 1;
		while ((j <= final) && bool)
			{
				j += ((j < final) && funcSort (vec[j], vec[j+1]));
				bool = funcSort (vec[k], vec[j]);
				if (bool)
					{ 
						HeapSortInts_Swap (vec+k, vec+j);
						k = j; j = (k << 1) + 1;
					}
			}
	}

void HeapSortInts (int *vec, int dim, Ints_Sort_func funcSort)
	{
		int i, ind, medio = dim / 2;

		for (i=1; i<=medio; i++)
			HeapSortInts_Adjust (vec, (medio-i), (dim-1), funcSort);
		for (i=0; i<(dim-1); i++)
			{
				ind = ((dim-2) - i);
				HeapSortInts_Swap (vec, vec+(ind+1));
				HeapSortInts_Adjust(vec, 0, ind, funcSort);
			}
	}

int FindSortInt (int num, int *vec, int dim, Ints_Sort_func funcSort)
	{
		int bool = 1, izq = 0, der = dim - 1, cnt = 0; 

		while ((bool) && (izq <= der))
			{
				cnt = (izq + der) / 2;
				if (num == vec[cnt])
					bool = 0;
				else if (funcSort (num, vec[cnt])) der = cnt - 1;
				else izq = cnt + 1;
			}
		if (bool) cnt = -1;
		else if (cnt > 0)
			{
				if (funcSort (0,1))
					{ while ((cnt >= 0) && (vec[cnt] == num)) cnt--; cnt++; }
				else
					{ while ((cnt < dim) && (vec[cnt] == num)) cnt++; cnt--; }
			}
		return cnt;
	}
	
int AddInts (int *vec, int num) 
	{
		int i, *pi = vec, aux = 0;
		for (i=0; i<num; i++) { aux += *pi; pi++; }
		return aux;
	}

int AddPermuteInts (int *vec, int *perm, int num) 
	{
		int i, *pi = perm, aux = 0;
		
		for (i=0; i<num; i++) { aux += vec[*pi]; pi++; }
		return aux;
	}

int PermuteInts (int *vec, int *perm, int num)
	{
		int i, *pi = vec;
		
		for (i=0; i<num; i++) { *pi = perm[*pi]; pi++; }
	}

int ComputeInvPermutation (int *perm, int *iperm, int num)
	{
		int i, *pi1 = perm;
		
		for (i=0; i<num; i++) { iperm[*(pi1++)] = i; }
	}

int PrintInts (int *vint, int num)
	{
		int i, *pi = vint;

		for (i=0; i<num; i++) printf ("%d ", *(pi++));
		printf ("\n");
	}

int FPrintInts (FILE * f, int *vint, int num)
	{
		int i, *pi = vint;

		for (i=0; i<num; i++) fprintf (f, "%d ", *(pi++));
		fprintf (f, "\n");
	}

int RemoveInts (int **vint)
	{ free (*vint); *vint = NULL; }

int CreateDoubles (double **vdouble, int num)
	{
		if ((*vdouble = (double *) malloc (sizeof(double)*num)) == NULL)
			{ printf ("Memory Error (CreateDoubles(%d))\n", num); exit (1); }
	}

int InitDoubles (double *vdouble, int n, double frst, double incr) 
	{
		int i; 
		double *p1 = vdouble, num = frst;

		for (i=0; i<n; i++) 
			{ *(p1++) = num; num += incr; }
	}
		
int CopyDoubles (double *src, double *dst, int n)
	{ 
		memmove (dst, src, sizeof(double) * n);
	}

void GetDoubleFromString (char *string, double *pdbl, int numC)
	{
		int j, k, exp, neg;
		double num, frac;
		char *pchar = string;

		j = 0; exp = 0; neg = 0; num = 0.0; frac = 1.0;
		while ((j < numC) && ((*pchar < '0') || (*pchar > '9')) && 
					 (*pchar != '+') && (*pchar != '-') && (*pchar != '.')) { j++; pchar++; }
		if (j < numC)
			{
				if ((*pchar == '+') || (*pchar == '-'))
					{ neg = (*pchar == '-'); j++; pchar++; }
				if (j < numC)
					{
						if (*pchar != '.')
							while ((j < numC) && (*pchar >= '0') && (*pchar <= '9'))
								{ num = num * 10 + (*pchar - 48); j++; pchar++; }
						if (j < numC)
							{
								if (*pchar == '.')
									{
										j++; pchar++; 
										while ((j < numC) && (*pchar >= '0') && (*pchar <= '9'))
											{ frac /= 10; num += (*pchar-48) * frac; j++; pchar++; }
									}
								if (neg) num = -num;
								if (j < numC)
									{
										if ((*pchar == 'e') || (*pchar == 'E') || (*pchar == 'd') || (*pchar == 'D'))
											{
												neg = 0; j++; pchar++; 
												if (j < numC)
													{
														if ((*pchar == '+') || (*pchar == '-'))
															{ neg = (*pchar == '-'); j++; pchar++; }
														if (j < numC)
															{
																while ((j < numC) && (*pchar >= '0') && 
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

void GetDoublesFromString (char *string, double *vec, int numN, int numC)
	{
		int i;
		double *paux = vec;
		char *pchar = string;

		for (i=0; i<numN; i++)
			{ GetDoubleFromString (pchar, (paux++), numC); pchar += numC; }
	}

void VvecDoubles (double alfa, double *src1, double *src2, double beta, double *dst, int num)
  {
    int i;

    for (i=0; i<num; i++)
      { *dst = (beta * *dst) + (alfa * *(src1++) * *(src2++)); dst++; }
  }

int PrintDoubles (double *vdouble, int num)
	{
		int i;
		double *pd = vdouble;

		for (i=0; i<num; i++) printf ("%e ", *(pd++));
		printf ("\n");
	}

int FPrintDoubles (FILE *f, double *vdouble, int num)
	{
		int i;
		double *pd = vdouble;

		for (i=0; i<num; i++) fprintf (f, "%20.14g ", *(pd++));
		fprintf (f, "\n");
	}

int RemoveDoubles (double **vdouble)
	{ free (*vdouble); *vdouble = NULL; }

