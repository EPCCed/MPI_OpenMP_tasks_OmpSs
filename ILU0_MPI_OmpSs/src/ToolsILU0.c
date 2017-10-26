#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "reloj.h"
#include <InputOutput.h>
#include "SparseSymmetricNew.h"
#include "ToolsILU0.h"

/********************************************************************************************/

void FactorILU0SparseMatrix (ptr_ILU0Factor pFact, paramFactor parFac,
                             int task, int leaf) {
	int i, j, k, i1, i2, j1, j2, cond, ierr = 0;
  int dimM, nnzM, indexM, dimL, nnzL, nnzF, nnzS;
  int *headL, *ptrs, *sizF;
  double d, val;
  SparseMatrix sprM, sprLDU, sprF, sprS;
 
   indexM = pFact->indM; headL = pFact->headL; dimL = headL[1];
   sprM = pFact->sprM; dimM = sprM.dim1; nnzM = sprM.vptr[dimM] - indexM;
 
   CreateInts (&ptrs, dimM); CopyInts (sprM.vptr, ptrs, dimM);
   CreateInts (&sizF, dimM); InitInts (sizF, dimM, 0, 0);
 
   if (leaf) {
   	RestrictUpperSparseMatrices (sprM, indexM, &sprLDU, indexM);
    HeapSortSparseMatrix (sprLDU, indexM, HeapSortSparseVector_IncrPos);
   } else {
	 	CreateSparseMatrix (&sprLDU, indexM, dimM, dimM, nnzM, 0);
    CopyInts    (sprM.vptr, sprLDU.vptr, dimM+1);
    CopyInts    (sprM.vpos, sprLDU.vpos, nnzM);
    CopyDoubles (sprM.vval, sprLDU.vval, nnzM);
   }
	 for (i=0; (i<dimL && ierr==0); i++) {
	 	if (sprLDU.vpos[sprLDU.vptr[i]] != (i+indexM)) {
      ierr = 1;  
    } else {
      cond = 1;
      d = sprLDU.vval[sprLDU.vptr[i]];
      for (j=sprLDU.vptr[i]+1; j<sprLDU.vptr[i+1]; j++) {
        if (cond && ((sprLDU.vpos[j]-indexM) >= dimL)) {
          ptrs[i] = j; cond = 0;
        }   
        sprLDU.vval[j] /= d;
      }   
      if (cond) ptrs[i] = sprLDU.vptr[i+1];
      sizF[i] = sprLDU.vptr[i+1] -  ptrs[i];
      d = sprLDU.vval[sprLDU.vptr[i]];
      for (j=sprLDU.vptr[i]+1; j<sprLDU.vptr[i+1]; j++) {
        k  = sprLDU.vpos[j]; val = sprLDU.vval[j  ] * d;
        i1 = sprLDU.vptr[i]; i2  = sprLDU.vptr[i+1];
        j1 = sprLDU.vptr[k]; j2  = sprLDU.vptr[k+1];
        while ((i1 < i2) && (j1 < j2)) {
					 if (sprLDU.vpos[i1] == sprLDU.vpos[j1]) {
					 		sprLDU.vval[j1] -= val * sprLDU.vval[i1];
            i1++; j1++;
          } else if (sprLDU.vpos[i1] < sprLDU.vpos[j1])
            i1++;
          else
            j1++;
        }
				     }
    }
  }
  if (ierr == 0) {
    nnzF = AddInts (sizF, dimL);
    nnzS = sprLDU.vptr[dimM] - sprLDU.vptr[dimL];
    nnzL = nnzM - (nnzF + nnzS);
	  // Creating rectangular matrix F
	   CreateSparseMatrix (&sprF, indexM, dimL, dimM-dimL, nnzF, 0);
    j = 0;
    for (i=0; i<dimL; i++) {
      sprF.vptr[i+1] = sprLDU.vptr[i+1] - ptrs[i];
      CopyShiftInts (sprLDU.vpos+ptrs[i], sprF.vpos+j, sprF.vptr[i+1], -dimL);
      CopyDoubles   (sprLDU.vval+ptrs[i], sprF.vval+j, sprF.vptr[i+1]);
      j += sprF.vptr[i+1];
    }
    *(sprF.vptr) = indexM; TransformLengthtoHeader (sprF.vptr, dimL);

    // Creating nonfactorized block S
     CreateSparseMatrix (&sprS, indexM, dimM-dimL, dimM-dimL, nnzS, 0);
    CopyShiftInts (sprLDU.vptr+dimL, sprS.vptr, dimM-dimL+1, -sprLDU.vptr[dimL]+indexM);
    CopyShiftInts (sprLDU.vpos+sprLDU.vptr[dimL]-indexM, sprS.vpos, nnzS, -dimL);
    CopyDoubles   (sprLDU.vval+sprLDU.vptr[dimL]-indexM, sprS.vval, nnzS);

    // Adjusting LDU
    InitInts (ptrs, dimL, 0, 1);  InitInts (ptrs+dimL, dimM-dimL, -1, 0);
    PermuteColsWithNegSparseMatrix (sprLDU, indexM, ptrs);
    sprLDU.dim1 = sprLDU.dim2 = dimL;
    sprLDU.vpos -= (dimM-dimL); CopyInts    (sprLDU.vpos+dimM-dimL, sprLDU.vpos, nnzL);
    ReallocSparseMatrix (&sprLDU);
	}
	pFact->ierr = ierr;
  pFact->sprLDU = sprLDU; 
  pFact->UF = sprF; pFact->sprF = sprS; pFact->indF = indexM;
  RemoveInts (&sizF); RemoveInts (&ptrs);
}

/******************************************************************************************************/

// This routine applies the preconditioner, going up.
// // void ILU0ResolutionUp (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf, int task) {
void ILU0ResolutionUp (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf) {
  // Definition of the local variables
	double *ptrI, *ptrO, *ptrA;
  double tle1, tle2, tlu1, tlu2;
  
  // Parameters validation
  if ((pFact == NULL) || (vidI < 0) || (vidO < 0) || (vidA < 0)) {
  	printf ("Incorrect parameters in ILU0ResolutionUp (%d,%d,%d)\n",
              vidI, vidO, vidA);
    PrintTrace (); exit (-1);
  }
  reloj (&tle1, &tlu1);
  // Definition of the vectors required in the computation
  ptrI = IdentifyVectorResolution (pFact, vidI);
  ptrO = IdentifyVectorResolution (pFact, vidO);
  ptrA = IdentifyVectorResolution (pFact, vidA);
  // Resolution of the local system
	{
    int i, j, indexM, dimL, dimUF;
    double aux;
    SparseMatrix sprLDU, UF;

    indexM = pFact->indM; 
		sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    CopyDoubles (ptrI, ptrA, dimL);
		sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    for (i=0; i<dimL; i++) {
      if (sprLDU.vpos[sprLDU.vptr[i]] == i+indexM) {
        ptrO[i] = ptrA[i];   
      } else {
        printf ("ILU0ResolutionUp Error in %d row\n", i);
      }
      for (j=sprLDU.vptr[i]+1; j<sprLDU.vptr[i+1]; j++) {
        ptrA[sprLDU.vpos[j]-indexM] -= ptrO[i] * sprLDU.vval[j];
      }
    }
	 	UF      = pFact->UF   ; dimUF = UF.dim2;
    CopyDoubles (ptrI+dimL, ptrO+dimL, dimUF);
    for (i=0; i<dimL; i++) {
      for (j=UF.vptr[i]; j<UF.vptr[i+1]; j++) {
        ptrO[dimL+UF.vpos[j]-indexM] -= ptrO[i] * UF.vval[j];
      }
    }
	}
	// Store the computation cost of the routine
	reloj (&tle2, &tlu2);
	pFact->tLoc[0][TLSOLU] += tle2-tle1; pFact->tLoc[1][TLSOLU] += tlu2-tlu1;
}

// This routine applies the preconditioner, in the top of the tree.
void ILU0ResolutionUpDown (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf) {
  // Definition of the local variables
  double *ptrI, *ptrO, *ptrA;
  double tle1, tle2, tlu1, tlu2;
  
  // Parameters validation
  if ((pFact == NULL) || (vidI < 0) || (vidO < 0) || (vidA < 0)) {
  	printf ("Incorrect parameters in ILU0ResolutionUpDown (%d,%d,%d)\n",
              vidI, vidO, vidA);
    PrintTrace (); exit (-1);
  }
  reloj (&tle1, &tlu1);
  // Definition of the vectors required in the computation
  ptrI = IdentifyVectorResolution (pFact, vidI);
  ptrO = IdentifyVectorResolution (pFact, vidO);
  ptrA = IdentifyVectorResolution (pFact, vidA);
  // Resolution of the local system
   {
    int i, j, indexM, dimL, dimUF;
    double aux;
    SparseMatrix sprLDU, UF;

    indexM = pFact->indM; 
		sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    CopyDoubles (ptrI, ptrA, dimL);

    // Computing the leading linear system L_B * y_B = r_B
    sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    for (i=0; i<dimL; i++) {
      if (sprLDU.vpos[sprLDU.vptr[i]] == i+indexM) {
        ptrA[i] = ptrA[i]; //  
      } else {
        printf ("ILU0ResolutionUpDown Error in %d row\n", i);
      }
      for (j=sprLDU.vptr[i]+1; j<sprLDU.vptr[i+1]; j++) {
        ptrA[sprLDU.vpos[j]-indexM] -= ptrA[i] * sprLDU.vval[j];
      }
    }

    // Scaling the input vector y_B = inv(D_B) * y_B
		sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    for (i=0; i<dimL; i++) {
      if (sprLDU.vpos[sprLDU.vptr[i]] == i+indexM) {
        ptrA[i] = ptrA[i] / sprLDU.vval[sprLDU.vptr[i]];
      } else {
        printf ("ILU0ResolutionUpDown Error in %d row\n", i);
      }
    }
		for (i=dimL-1; i>=0; i--) {
      aux = 0;
      for (j=sprLDU.vptr[i]+1; j<sprLDU.vptr[i+1]; j++) {
        aux += ptrO[sprLDU.vpos[j]-indexM] * sprLDU.vval[j];
      }
      if (sprLDU.vpos[sprLDU.vptr[i]] == i+indexM) {
        ptrO[i] = (ptrA[i]-aux); 
      } else {
        printf ("ILU0ResolutionUpDown Error in %d row\n", i);
      }
    }
	}
	// Store the computation cost of the routine
  reloj (&tle2, &tlu2);
  pFact->tLoc[0][TLSOLU] += tle2-tle1; pFact->tLoc[1][TLSOLU] += tlu2-tlu1;
}

// This routine applies the preconditioner, going down.
void ILU0ResolutionDown (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf) {
  // Definition of the local variables
  double *ptrI, *ptrO, *ptrA;
  double tle1, tle2, tlu1, tlu2;
  
  // Parameters validation
  if ((pFact == NULL) || (vidI < 0) || (vidO < 0) || (vidA < 0)) {
  	printf ("Incorrect parameters in ILU0ResolutionDown (%d,%d,%d)\n",
              vidI, vidO, vidA);
    PrintTrace (); exit (-1);
  }
  reloj (&tle1, &tlu1);
  // Definition of the vectors required in the computation
  ptrI = IdentifyVectorResolution (pFact, vidI);
  ptrO = IdentifyVectorResolution (pFact, vidO);
  ptrA = IdentifyVectorResolution (pFact, vidA);
  // Resolution of the local system
  {
    int i, j, indexM, dimL, dimUF;
    double aux;
    SparseMatrix sprLDU, UF;

    indexM = pFact->indM; 
		sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    CopyDoubles (ptrI, ptrA, dimL);
		// Scaling the input vector y_B = inv(D_B) * y_B
		sprLDU = pFact->sprLDU; dimL = sprLDU.dim1;
    for (i=0; i<dimL; i++) {
      if (sprLDU.vpos[sprLDU.vptr[i]] == i+indexM) {
        ptrA[i] = ptrA[i] / sprLDU.vval[sprLDU.vptr[i]];
      } else {
        printf ("ILU0ResolutionDown Error in %d row\n", i);
      }
    }
	 	// Copy the processed block zc = yc
	  UF      = pFact->UF   ; dimUF = UF.dim2;
	  CopyDoubles (ptrI+dimL, ptrO+dimL, dimUF);
  	// Updating the next level vector y_B = y_B - U_F * z_C 
  	for (i=0; i<dimL; i++) {
      aux = 0;
      for (j=UF.vptr[i]; j<UF.vptr[i+1]; j++) {
        aux += ptrO[dimL+UF.vpos[j]-indexM] * UF.vval[j];
      }
      ptrA[i] = ptrA[i] - aux;
    }
    // Computing the leading linear system U_B * z_B = y_B
     for (i=dimL-1; i>=0; i--) {
      aux = 0;
      for (j=sprLDU.vptr[i]+1; j<sprLDU.vptr[i+1]; j++) {
        aux += ptrO[sprLDU.vpos[j]-indexM] * sprLDU.vval[j];
      }
      if (sprLDU.vpos[sprLDU.vptr[i]] == i+indexM) {
        ptrO[i] = (ptrA[i]-aux); 
      } else {
        printf ("ILU0ResolutionDown Error in %d row\n", i);
      }
    }
	}
	// Store the computation cost of the routine
  reloj (&tle2, &tlu2);
  pFact->tLoc[0][TLSOLD] += tle2-tle1; pFact->tLoc[1][TLSOLD] += tlu2-tlu1;
}
