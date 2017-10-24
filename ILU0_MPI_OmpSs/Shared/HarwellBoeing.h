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
#ifndef HarwellBoeing

#define HarwellBoeing 1

extern void ReadHeaderFileHB (char *nameFile, int *numF, int *numC, int *numE);

extern void ReadPtrHB (char *nameFile, int *vptr, int FtoC);

#endif
