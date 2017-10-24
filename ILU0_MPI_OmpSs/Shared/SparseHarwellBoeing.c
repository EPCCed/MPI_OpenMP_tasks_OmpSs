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
#include "SparseMatrices.h"
#include "HarwellBoeing.h"
#include "SparseHarwellBoeing.h"

#define length1 82
#define length2 82

void ReadSparseIndexValuesHB (char *nameFile, int *vpos, double *vval, int FtoC)	
	{
		FILE *file;
		char string[length1]; 
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int lines[5], dim[4], formats[10];

		file = OpenFile (nameFile, "r");

		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0);

		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		ReadStringFile (file, string, length1);
		GetFormatsFromString (string, formats, 2, 16);
		GetFormatsFromString ((string+32), (formats+4), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) ReadStringFile (file, string, length1);

		for (i = 0; i < lines[1]; i++)
			ReadStringFile (file, string, length2);

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2);
				k = (dim[2] - j);
				if (k > formats[2]) k = formats[2];
				GetIntsFromString (string, (vpos+j), k, formats[3], shft);
				j+=formats[2];
			}

		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2);
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetDoublesFromString (string, (vval+j), k, formats[5]);
				j+=formats[4];
			}

		fclose (file);
	}

void ReadSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC)
	{
		FILE *file;
		char string[length1]; 
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int *vptr = spr->vptr, *vpos = spr->vpos;
		double *vval = spr->vval; 
		int lines[5], dim[4], formats[10];

		file = OpenFile (nameFile, "r");

		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0);

		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		ReadStringFile (file, string, length1);
		GetFormatsFromString (string, formats, 2, 16);
		GetFormatsFromString ((string+32), (formats+4), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) ReadStringFile (file, string, length1);

		j = 0;
		for (i = 0; i < lines[1]; i++)
			{
				ReadStringFile (file, string, length2);
				k = ((dim[0] + 1) - j);
				if (k > formats[0]) k = formats[0];
				GetIntsFromString (string, (vptr+j), k, formats[1], shft);
				j+=formats[0];
			}

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2);
				k = (dim[2] - j);
				if (k > formats[2]) k = formats[2];
				GetIntsFromString (string, (vpos+j), k, formats[3], shft);
				j+=formats[2];
			}

		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2);
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetDoublesFromString (string, (vval+j), k, formats[5]);
				j+=formats[4];
			}

		fclose (file);
	}

void CreateSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC)
	{
		FILE *file;
		char string[length1]; 
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int *vptr = NULL, *vpos = NULL;
		double *vval = NULL; 
		int lines[5], dim[4], formats[10];

		file = OpenFile (nameFile, "r");
		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0); 
		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		CreateSparseMatrix (spr, dim[0], dim[1], dim[2], 0);
		vptr = spr->vptr; vpos = spr->vpos; vval = spr->vval; 

		ReadStringFile (file, string, length1);
		GetFormatsFromString (string, formats, 2, 16);
		GetFormatsFromString ((string+32), (formats+4), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) ReadStringFile (file, string, length1);

		j = 0;
		for (i = 0; i < lines[1]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = ((dim[0] + 1) - j);
				if (k > formats[0]) k = formats[0];
				GetIntsFromString (string, (vptr+j), k, formats[1], shft);
				j+=formats[0];
			}

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[2]) k = formats[2];
				GetIntsFromString (string, (vpos+j), k, formats[3], shft);
				j+=formats[2];
			}

		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetDoublesFromString (string, (vval+j), k, formats[5]);
				j+=formats[4];
			}

		fclose (file);
	}

void CreateSparseMatrixHB2 (char *nameFile, ptr_SparseMatrix spr, int FtoC)
	{
		FILE *file;
		char string[length1]; 
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int *vptr = NULL, *vpos = NULL;
		double *vval = NULL; 
		int lines[5], dim[4], formats[12];

		file = OpenFile (nameFile, "r");

		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0); 
		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		CreateSparseMatrix (spr, dim[0], dim[1], dim[2], 0);
		vptr = spr->vptr; vpos = spr->vpos; vval = spr->vval; 

		ReadStringFile (file, string, length1);
		GetFormatsFromString2 (string, formats, 2, 16);
		GetFormatsFromString2 ((string+32), (formats+6), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) ReadStringFile (file, string, length1);

		j = 0;
		for (i = 0; i < lines[1]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = ((dim[0] + 1) - j);
				if (k > formats[1]) k = formats[1];
				GetIntsFromString (string+formats[0], (vptr+j), k, formats[2], shft);
				j+=formats[1];
			}

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetIntsFromString (string+formats[3], (vpos+j), k, formats[5], shft);
				j+=formats[4];
			}

		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[6]) k = formats[6];
				GetDoublesFromString (string, (vval+j), k, formats[7]);
				j+=formats[6];
			}

		fclose (file);
	}

void CreateSparseMatrixHB3 (char *nameFile, ptr_SparseMatrix spr, double **b, double **x, int FtoC)
	{
		FILE *file;
		char string[length1]; 
		char c0 = ' ', c1 = ' ', c2 = ' ';
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int *vptr = NULL, *vpos = NULL;
		double *vval = NULL; 
		int lines[5], dim[6], formats[12];

		file = OpenFile (nameFile, "r");

		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0); 
		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		CreateSparseMatrix (spr, dim[0], dim[1], dim[2], 0);
		vptr = spr->vptr; vpos = spr->vpos; vval = spr->vval; 

		ReadStringFile (file, string, length1);
		GetFormatsFromString2 (string, formats, 2, 16);
		GetFormatsFromString2 ((string+32), (formats+6), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) 
			{
				ReadStringFile (file, string, length1);
				c0 = string[0]; c1 = string[1]; c2 = string[2];
				GetIntsFromString ((string+14), dim+4, 2, 14, 0);
			}

		PrintInts (lines, 5); PrintInts (dim, 4+(2*(lines[4] > 0))); PrintInts (formats, 12);

		j = 0;
		for (i = 0; i < lines[1]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = ((dim[0] + 1) - j);
				if (k > formats[1]) k = formats[1];
				GetIntsFromString (string+formats[0], (vptr+j), k, formats[2], shft);
				j+=formats[1];
			}

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetIntsFromString (string+formats[3], (vpos+j), k, formats[5], shft);
				j+=formats[4];
			}

		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[6]) k = formats[6];
				GetDoublesFromString (string, (vval+j), k, formats[7]);
				j+=formats[6];
			}

		if ((lines[4] > 0) && ((c0 == 'f') || (c0 == 'F')))
			{
				CreateDoubles (b, dim[0]*dim[4]); 
				if ((c2 == 'x') || (c2 == 'X')) CreateDoubles (x, dim[0]); 
				else	*x = NULL;
				j = 0; vval = *b;
				for (i = 0; i < lines[4]; i++)
					{
						ReadStringFile (file, string, length2); 
						k = (dim[0] - j);
						if (k > formats[9]) k = formats[9];
						GetDoublesFromString (string, (vval+j), k, formats[10]);
						j+=formats[9]; if (j >= dim[0]*dim[4]) { vval = *x; j = 0; }
					}
			}
		else 
			{ *x = NULL; *b = NULL; }

		fclose (file);
	}

void WriteSparseMatrixHB (char *nameFile, SparseMatrix spr, int numI, int numR, 
													char *tag, int CtoF) 
	{
		FILE *file;
		int i, j, numE, M = spr.dim1, N = spr.dim2, shft = (CtoF)?1:0, lines[4];
		char formI[8], formR[12];
		int *vptr = spr.vptr, *vpos = spr.vpos;
		double *vval = spr.vval; 

		if ((spr.dim1 > 0) && (spr.dim2 > 0)) {
			file = OpenFile (nameFile, "w");

			PrintSpaceFile (file, 72);
			fprintf (file, "%s\n", tag);

			numE = vptr[M];
			lines[1] = ((M + 1)	/ numI) + (((M + 1)	% numI) > 0);
			lines[2] = (numE / numI) + ((numE % numI) > 0);
			lines[3] = (numE / numR) + ((numE % numR) > 0);
			lines[0] = lines[1] + lines[2] + lines[3];
			fprintf (file, "%14d%14d%14d%14d%14d\n", 
			lines[0], lines[1], lines[2], lines[3], 0);

			fprintf (file, "RSA"); PrintSpaceFile (file, 11);
			fprintf (file, "%14d%14d%14d%14d\n", N, M, numE, 0);

			fprintf (file, "(%dI%d)", numI, 80 / numI); PrintSpaceFile (file, 10);
			fprintf (file, "(%dI%d)", numI, 80 / numI); PrintSpaceFile (file, 10);
			fprintf (file, "(%dE%d.%d)", numR, (80 / numR), 40 / numR); 
			PrintSpaceFile (file, 11); fprintf (file, "(*)\n");
			sprintf (formI, "%%%dd", 80 / numI);
			i = (80/numR); j = i-7;
			sprintf (formR, "%%%d.%de", i, j);

			for (i = 0; i <= M; i++)
				{
					fprintf (file, formI, *(vptr++)+shft);
					if ((i % numI) == (numI-1)) fprintf (file, "\n");
				}
			if ((M % numI) != (numI-1)) fprintf (file, "\n");

			for (i = 0; i < numE; i++)
				{
					fprintf (file, formI, *(vpos++)+shft);
					if ((i % numI) == (numI-1)) fprintf (file, "\n");
				}
			if ((numE % numI) != 0) fprintf (file, "\n");

			for (i = 0; i < numE; i++)
				{
					fprintf (file, formR, *(vval++));
					if ((i % numR) == (numR - 1)) fprintf (file, "\n");
				}
			if ((numE % numR) != 0) fprintf (file, "\n");

			fclose (file);
		}
	}

