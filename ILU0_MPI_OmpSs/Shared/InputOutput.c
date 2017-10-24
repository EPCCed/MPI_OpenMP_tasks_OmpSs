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
#include <execinfo.h>
#include "InputOutput.h"

#define screen 29

void ClearScreen ()
	{
		int i;
		for (i=0; i<screen; i++) printf ("\n");
	}

char ReadCharacter ()
	{
		char c;
		do 
			{ c = getchar (); } 
		while (((c < 'a') || (c > 'z')) && ((c < 'A') || (c > 'Z')));
		return c;
	}

char ReadAlfanum ()
	{
		char c;
		do 
			{ c = getchar (); } 
		while (((c < 'a') || (c > 'z')) && ((c < 'A') || (c > 'Z')) &&
					 ((c < '0') || (c > '9')));
		return c;
	}

unsigned int ReadUnsignedInt ()
	{
		char c;
		unsigned int num = 0;
		do 
			{ c = getchar (); } 
		while ((c < '0') || (c > '9'));
		do 
			{ num = num * 10 + (c - 48); c = getchar (); }
		while ((c >= '0') && (c <= '9'));
		return num;
	}

int ReadInt ()
	{
		char c;
		int num = 0, neg = 0;
		do 
			{ c = getchar (); } 
		while (((c < '0') || (c > '9')) && (c != '+') && (c != '-'));
		if ((c == '+') || (c == '-'))
			{ neg = (c == '-'); c = getchar (); }
		while ((c >= '0') && (c <= '9'))
			{ num = num * 10 + (c - 48); c = getchar (); }
		if (neg) num = -num;
		return num;
	}

float ReadUnsignedReal ()
	{
		char c;
		int i, exp = 0, neg = 0;
		float num = 0, frac = 1;

		do 
			{ c = getchar (); } 
		while (((c < '0') || (c > '9')) && (c != '.'));
		if (c != '.')
			while ((c >= '0') && (c <= '9'))
				{ num = num * 10 + (c - 48); c = getchar (); }
		if (c == '.')
			{ 
				c = getchar ();
				while ((c >= '0') && (c <= '9'))
					{ frac /= 10; num += (c - 48) * frac; c = getchar (); }
			}
		if ((c == 'e') || (c == 'E') || (c == 'd') || (c == 'D'))
			{
				neg = 0; c = getchar ();
				if ((c == '+') || (c == '-'))
					{ neg = (c == '-'); c = getchar (); }
				while ((c >= '0') && (c <= '9'))
					{ exp = exp * 10 + (c - 48); c = getchar (); }
				if (neg) exp = -exp;
				for (i=0; i<exp; i++) num *= 10;
				for (i=0; i>exp; i--) num /= 10;
			}
		return num;
	}

float ReadReal ()
	{
		char c;
		int i, exp = 0, neg = 0;
		float num = 0, frac = 1;
		do 
			{ c = getchar (); } 
		while (((c < '0') || (c > '9')) && (c != '+') && (c != '-') && (c != '.'));
		if ((c == '+') || (c == '-'))
			{ neg = (c == '-'); c = getchar (); }
		if (c != '.')
			while ((c >= '0') && (c <= '9'))
				{ num = num * 10 + (c - 48); c = getchar (); }
		if (c == '.')
			{ 
				c = getchar ();
				while ((c >= '0') && (c <= '9'))
					{ frac /= 10; num += (c - 48) * frac; c = getchar (); }
			}
		if (neg) num = -num;
		if ((c == 'e') || (c == 'E') || (c == 'd') || (c == 'D'))
			{
				neg = 0; c = getchar ();
				if ((c == '+') || (c == '-'))
					{ neg = (c == '-'); c = getchar (); }
				while ((c >= '0') && (c <= '9'))
					{ exp = exp * 10 + (c - 48); c = getchar (); }
				if (neg) exp = -exp;
				for (i=0; i<exp; i++) num *= 10;
				for (i=0; i>exp; i--) num /= 10;
			}
		return num;
	}

double ReadUnsignedDouble ()
	{
		char c;
		int i, exp = 0, neg = 0;
		double num = 0, frac = 1;
		do 
			{ c = getchar (); } 
		while (((c < '0') || (c > '9')) && (c != '.'));
		if (c != '.')
			while ((c >= '0') && (c <= '9'))
				{ num = num * 10 + (c - 48); c = getchar (); }
		if (c == '.')
			{ 
				c = getchar ();
				while ((c >= '0') && (c <= '9'))
					{ frac /= 10; num += (c - 48) * frac; c = getchar (); }
			}
		if ((c == 'e') || (c == 'E') || (c == 'd') || (c == 'D'))
			{
				neg = 0; c = getchar ();
				if ((c == '+') || (c == '-'))
					{ neg = (c == '-'); c = getchar (); }
				while ((c >= '0') && (c <= '9'))
					{ exp = exp * 10 + (c - 48); c = getchar (); }
				if (neg) exp = -exp;
				for (i=0; i<exp; i++) num *= 10;
				for (i=0; i>exp; i--) num /= 10;
			}
		return num;
	}

double ReadDouble ()
	{
		char c;
		int i, exp = 0, neg = 0;
		double num = 0, frac = 1;
		do 
			{ c = getchar (); } 
		while (((c < '0') || (c > '9')) && (c != '+') && (c != '-') && (c != '.'));
		if ((c == '+') || (c == '-'))
			{ neg = (c == '-'); c = getchar (); }
		if (c != '.')
			while ((c >= '0') && (c <= '9'))
				{ num = num * 10 + (c - 48); c = getchar (); }
		if (c == '.')
			{ 
				c = getchar ();
				while ((c >= '0') && (c <= '9'))
					{ frac /= 10; num += (c - 48) * frac; c = getchar (); }
			}
		if (neg) num = -num;
		if ((c == 'e') || (c == 'E') || (c == 'd') || (c == 'D'))
			{
				neg = 0; c = getchar ();
				if ((c == '+') || (c == '-'))
					{ neg = (c == '-'); c = getchar (); }
				while ((c >= '0') && (c <= '9'))
					{ exp = exp * 10 + (c - 48); c = getchar (); }
				if (neg) exp = -exp;
				for (i=0; i<exp; i++) num *= 10;
				for (i=0; i>exp; i--) num /= 10;
			}
		return num;
	}

void WaitChar ()
	{ char c; fflush (stdin); printf ("\nPress a key --> "); c = ReadAlfanum (); printf ("\n"); }

int YesOrNot (char *str)
	{
		char c; 
		do
			{ printf ("\n%s (Y/N) --> ", str); c = ReadCharacter (); }
		while ((c != 'y') && (c != 'Y') && (c != 'n') && (c != 'N'));
		printf ("\n");
		return ((c == 'y') || (c == 'Y'));
	}

int SelectOption (int min, int max, char *str)
	{
		int i;
		do 
			{ printf ("\n%s --> ", str); i = ReadInt (); }
		while ((i < min) || (i > max));
		printf ("\n");
		return i;
	}

void ChangeBinary (void *c, int size)
	{
		int i;
		char b, *c1 = c, *c2 = (char *) c + (size - 1);
		for (i=0; i<size; i+=2)
			{ b = *c1; *c1 = *c2; *c2 = b; c1++; c2 -= 2; }
	}

void ReadBinaryFile (FILE *f, void *c, int size)
  {
		if (fread (c, size, 1, f) == 0) 
			{ printf("Read Error %d in ReadBinaryFile\n",size); exit(1); }
#ifndef DoChanges
		ChangeBinary (c, size);
#endif
	}

void WriteBinaryFile (FILE *f, void *c, int size)
  {
#ifndef DoChanges
		ChangeBinary (c, size);
#endif
		if (fwrite(c, size, 1, f) == 0) 
			{ printf("Write Error %d in WriteBinaryFile\n",size); exit(1); }
#ifndef DoChanges
		ChangeBinary (c, size);
#endif
	}

FILE *OpenFile (char *name, char *attr)
	{
    FILE *fich;
//    printf ("Opening file %s\n", name);
    if ((fich = fopen (name, attr)) == NULL)
      { printf ("File %s not exists \n", name); exit(1); }
		return fich;
	}

void ReadStringFile (FILE *file, char *string, int length)
	{
		char *s = NULL;
		if ((s = fgets (string, length, file)) == NULL)
			{ printf ("Read Error in ReadStringFile\n"); exit (1); }
	}

void PrintSpaceFile (FILE *file, int num)
	{
		int i;
		for (i=0; i<num; i++) fprintf (file, " ");
	}

void PrintTrace (void) {
	void *array[100];
	char **strings;
	int size, i;
     
	size = backtrace (array, 100);
	strings = backtrace_symbols (array, size);
 
	if (strings == NULL) {	
		; 
	} else {
		printf ("Obtained %d stack frames.\n", size);
		for (i = 0; i < size; i++)
			printf ("%s\n", strings[i]);
		free (strings);
	}
}

