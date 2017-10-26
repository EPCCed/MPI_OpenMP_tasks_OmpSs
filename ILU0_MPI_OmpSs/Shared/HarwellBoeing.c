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
#include "HarwellBoeing.h"

#define length1 82
#define length2 82

void ReadHeaderFileHB (char *nameFile, int *numF, int *numC, int *numE) 
	{
		FILE *file;
		char string[length1];
		int num[4];

		file = OpenFile (nameFile, "r");

		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		fclose (file);
		GetIntsFromString ((string+14), num, 4, 14, 0);
		*numF	= num[0]; *numC	= num[1]; *numE = num[2];
	}

void ReadPtrHB (char *nameFile, int *vptr, int FtoC) 
	{
		FILE *file;
		char string[length1]; 
		int lines[5], dim[4], formats[10];
		int i, j, k = 0, shft = (FtoC)?-1:0;

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

		fclose (file);
	}

