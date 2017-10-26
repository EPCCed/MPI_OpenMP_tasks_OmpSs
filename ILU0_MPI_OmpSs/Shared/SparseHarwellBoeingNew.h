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
#ifndef SparseHarwellBoeingNew

#define SparseHarwellBoeingNew 1

#include "SparseMatricesNew.h"

extern void ReadSparseIndexValuesHB (char *nameFile, int *vpos, double *vval, int FtoC);

extern void ReadSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC);

extern void CreateSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC);

extern void CreateSparseMatrixHB2 (char *nameFile, ptr_SparseMatrix spr, int FtoC);

extern void WriteSparseMatrixHB (char *nameFile, SparseMatrix spr, int numI, int numR, 
													char *tag, int CtoF);

#endif
