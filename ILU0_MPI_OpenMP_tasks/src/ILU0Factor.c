#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "InputOutput.h"
#include "ILU0Factor.h"
#include "SparseHarwellBoeingNew.h"

/*********************************************************************************/

// This routine verify the number of parameters and its contents
void VerificationParameters (int argc, char **argv, int *nleaves, paramFactor *parFac) {
	int NPARTS, NTHREADS;
	double ELBOW, CONDEST, DROPTOL;
     
	// First, the number of paramters is checked
	if (argc != NUMBER_PARAMS) {
		printf("usage '%s <drop tol.> <bound for L^{-1},U^{-1}> <elbow space> <nthreads> <nparts> <matrix>'\n",argv[0]);
		exit(0);
	}
	// Getting the data from the command line
	NPARTS   = atoi(argv[argc-2]);
 	NTHREADS = atoi(argv[argc-3]);
 	ELBOW    = atof(argv[argc-4]);
 	CONDEST  = atof(argv[argc-5]);
 	DROPTOL  = atof(argv[argc-6]);
	// Verification of the values
	if ((NPARTS < 1) || (NTHREADS < 1) || (ELBOW <= 0) || (CONDEST <= 0) || (DROPTOL <= 0)) {
		printf("usage '%s <drop tol.> <bound for L^{-1},U^{-1}> <elbow space> <nthreads> <nparts> <matrix>'\n",argv[0]);
		printf("  where (droptol>0) && (condest>0) && (elbow>0) && (nparts>0) && (nparts>0)\n"); 
		exit(0);
	}

	// Definition of the number of leaves to be used in the computations
	*nleaves = NPARTS;
   
	// Definition of the parameter which are used on the computations
	parFac->nparts    = NPARTS      ;
	parFac->nthreads  = NTHREADS    ;
	parFac->maxit     = 500        ;
}


/*********************************************************************************/

// Create a vector of preconditioners with dim elements
ptr_ILU0Factor CreateILU0FactorVector (int dim) {
	// Definition of the local vectors and variables
	int i;
	SparseMatrix spr = {0 , 0, NULL, NULL, NULL};
	ptr_ILU0Factor vFact = NULL;

	// Parameters validation
	if (dim < 1) {
		printf ("Incorrect parameter in CreateILU0FactorVector (%d)\n", dim); 
		PrintTrace (); exit (-1);
	}
	// Create the vector of preconditioners
	if ((vFact = (ptr_ILU0Factor) malloc (dim * sizeof(ILU0Factor))) == NULL) {
      printf ("Memory Error (CreateILU0FactorVector (%d))\n", dim); exit (1); 
	}
  // Initialize the ILU0 structures
	vFact->UF    = vFact->sprLDU = spr;
	vFact->ierr  = 0; 
	// Initialize the vector of preconditioners
	vFact->indP  = vFact->indM     = vFact->indF  = vFact->indR = -1;
	vFact->perm  = NULL;
	vFact->mTab  = NULL;
	vFact->dimL  = vFact->dimT     = 0;
	vFact->dimX  = vFact->dimF     = 0;
	vFact->nlev  = 0; vFact->headL = NULL;  vFact->divs = NULL;
	vFact->sprM  = vFact->sprF     = spr ;
	vFact->dimM  = 0;
	vFact->mDia  = NULL;
	vFact->permM = NULL;
	vFact->sprR  = spr ;
	vFact->dimV  = 0;
	vFact->tskV  = -1;
	vFact->mPCG  = vFact->mVcL     = NULL ;
//	vFact->goup  = vFact->godown   = NULL;
	vFact->ntsks = NULL;
	vFact->tGlb  = vFact->tLoc     = NULL;
	for (i=1; i<dim; i++) vFact[i] = *vFact;
	// Return the created structure
	return vFact;
}

// Clean a vector of preconditioners, maintaining only the preconditioner,
// the permutation and the METIS structures.
void CleanILU0FactorVector (ptr_ILU0Factor *ppFact) {
	// Definition of the local vectors and variables
	int i;

	if (*ppFact != NULL) {
		// Free the components of the vector of preconditioner
		for (i=0; i<(*ppFact)[0].dimL; i++) {
			// Free the non-ILUPACK data
			RemoveSparseMatrix (&(*ppFact)[i].sprM); 
			RemoveMatrixDoubles (&(*ppFact)[i].mDia);
			RemoveSparseMatrix (&(*ppFact)[i].sprR); 
			if ((*ppFact)[i].tskV == i) RemoveMatrixDoubles (&(*ppFact)[i].mPCG); 
			else if ((*ppFact)[i].tskV >= 0) RemoveVectorPtrDoubles (&(*ppFact)[i].mPCG);
			if ((*ppFact)[i].tskV == i) RemoveMatrixDoubles (&(*ppFact)[i].mVcL); 
			else if ((*ppFact)[i].tskV >= 0) RemoveVectorPtrDoubles (&(*ppFact)[i].mVcL);
		}
	}
}
// Remove a vector of preconditioners 
void RemoveILU0FactorVector (ptr_ILU0Factor *ppFact) {
	// Definition of the local vectors and variables
	int i;

	if (*ppFact != NULL) {
		// Free the Metis structures
		RemoveMatrixInts (&(*ppFact)[0].mTab); RemoveInts (&(*ppFact)[0].perm); 
		RemoveMatrixDoubles (&(*ppFact)[0].tGlb);
		// Free the components of the vector of preconditioner
		for (i=0; i<(*ppFact)[0].dimL; i++) {
			// Free the non-ILU0 data
			if ((*ppFact)[i].sprM.vptr != (*ppFact)[i].sprR.vptr)
				RemoveSparseMatrix (&(*ppFact)[i].sprM); 
			RemoveMatrixDoubles (&(*ppFact)[i].mDia);
			RemoveDoubles (&(*ppFact)[i].divs); 
			RemoveInts (&(*ppFact)[i].permM); (*ppFact)[i].headL = NULL;
			RemoveSparseMatrix (&(*ppFact)[i].sprR); 
			if ((*ppFact)[i].tskV == i) RemoveMatrixDoubles (&(*ppFact)[i].mPCG); 
			else if ((*ppFact)[i].tskV >= 0) RemoveVectorPtrDoubles (&(*ppFact)[i].mPCG);
			if ((*ppFact)[i].tskV == i) RemoveMatrixDoubles (&(*ppFact)[i].mVcL); 
			else if ((*ppFact)[i].tskV >= 0) RemoveVectorPtrDoubles (&(*ppFact)[i].mVcL);
			// Free the ILU0 data
			if ((*ppFact)[i].sprLDU.vptr != NULL) RemoveSparseMatrix (&(*ppFact)[i].sprLDU);
			if ((*ppFact)[i].UF.vptr != NULL) RemoveSparseMatrix (&(*ppFact)[i].UF);

		}
		// Free the vector of preconditioners
		free (*ppFact); *ppFact = NULL;
	}
}

// Write a vector of preconditioners 
void WriteILU0FactorVector (ptr_ILU0Factor vFact, char *nameMatr, char *nameFact, char *nameScal) {
	// Definition of the global vectors and variables
	char nameFile[80], tag[4] = "JIA";
	int i, numI = 10, numR = 4;

	// Parameters validation
	if (vFact == NULL) {
		printf ("Incorrect parameters in WriteILU0FactorVector\n"); 
		PrintTrace (); exit (-1);
	}
	// Print the components of the vector of preconditioner
	for (i=0; i<vFact[0].dimL; i++) {
 		// Print the non-ILU0 data
		if (nameMatr != NULL) {
 			sprintf (nameFile, "%s_%-3.3d", nameMatr, i); // printf ("%s\n", nameFile);
 			WriteSparseMatrixHB (nameFile, vFact[i].sprM, numI, numR, tag, (1-vFact[i].indM));
		}
		if (nameFact != NULL) {
 			sprintf (nameFile, "%s_%-3.3d", nameFact, i); // printf ("%s\n", nameFile);
 			WriteSparseMatrixHB (nameFile, vFact[i].sprF, numI, numR, tag, (1-vFact[i].indF));
		}
		if (nameScal != NULL) {
 			sprintf (nameFile, "%s_%-3.3d", nameScal, i); // printf ("%s\n", nameFile);
		}
	}
}

// This routine returns the accumulation of the sizes of the leaves: 
// - The addition of the number of rows
// - The addition of the number of nonzeros
void GetSizesLeaves (ptr_ILU0Factor vFact, int *dimML, int *nnzML) {
	int i, dim = 0, nnz = 0;

	// Parameters validation
	if ((vFact == NULL) || (dimML == NULL) || (nnzML == NULL)) {
		printf ("Incorrect parameters in GetSizesLeavesILU0Factor\n"); 
		PrintTrace (); exit (-1);
	}
	// Print the components of the vector of preconditioner
	for (i=0; i<vFact[0].dimL; i++) {
		if (vFact->mTab[CHLD][i] == -1) {
			dim += vFact[i].dimM;
			if (vFact[i].sprM.dim1 > 0)
				nnz += vFact[i].sprM.vptr[vFact[i].dimM];
		}
	}
	*dimML = dim; *nnzML = nnz;
}

// This routine returns the accumulation of the sizes of the leaves: 
// - The addition of the number of rows
// - The addition of the number of nonzeros
void GetSizesLeavesPCG (ptr_ILU0Factor vFact, int *dimRL, int *nnzRL) {
	int i, dim = 0, nnz = 0;

	// Parameters validation
	if ((vFact == NULL) || (dimRL == NULL) || (nnzRL == NULL)) {
		printf ("Incorrect parameters in GetSizesLeavesILU0Factor\n"); 
		PrintTrace (); exit (-1);
	}
	// Print the components of the vector of preconditioner
	for (i=0; i<vFact[0].dimL; i++) {
		if (vFact->mTab[CHLD][i] == -1) {
			dim += vFact[i].dimM;
			if (vFact[i].sprR.dim1 > 0)
				nnz += vFact[i].sprR.vptr[vFact[i].dimM];
		}
	}
	*dimRL = dim; *nnzRL = nnz;
}

/*********************************************************************************/
// This routine returns the address of the vector vid related to the structure pFact
double *IdentifyVectorResolution (ptr_ILU0Factor pFact, int vid) {
	// Definition of the local vectors and variables
	double *ptr = NULL;
	int i = vid / 10, j = vid % 10;

	// Parameters validation
	if ((pFact == NULL) || (vid < 0)) {
		printf ("Incorrect parameters in IdentifyVectorResolution (%d)\n", vid); 
		PrintTrace (); exit (-1);
	}
	// The first vectors are related to the 
	switch (i) {
		case  0: ptr = ((j < SIZE_MAT_PCG)? pFact->mPCG[j]: NULL); break;
		case  1: ptr = ((j < SIZE_MAT_VCL)? pFact->mVcL[j]: NULL); break;
		default: break;
	}
	if (ptr == NULL) {
		printf ("Incorrect parameter in IdentifyVectorResolution (%d)\n", vid); 
		PrintTrace (); exit (-1);
	}
	// Return the pointer
	return ptr;
}

/****************************************************************************/

// This routine creates the vector tRes, whose size is equal the the number of leaves, 
// and copies the contents of a local timer to this new vector.
void GetLocalTimer (ptr_ILU0Factor vFact, int timer, double **tRes) {
	int i;
	double *times = NULL;

	if ((vFact == NULL) || (timer < 0) || (timer >= SIZE_TIM_LOC)) {
		printf ("Incorrect parameters in GetLocalTimer (%d)\n", timer);
		PrintTrace (); exit (-1);
	}
	CreateDoubles (&times, vFact->dimL);
	for (i=0; i<vFact->dimL; i++)
		if (vFact[i].tLoc == NULL)
      times[i] = 0.0;
    else
			times[i] = vFact[i].tLoc[0][timer];

	*tRes = times;
}

// This routine prints the information included in the counters delimited by
// nTim1 and nTim2, acording to the value of typeP:
// * < 1, it prints all the information, from the nodes to the final summarize.
// * < 2, it prints the information, from the threads to the final summarize.
// * < 3, it only prints the summarize.
// * < 4, it only computes the parameters.
// First, the information is accumulated by each level of the threads on the 
// matrix mSumThr, and after the addition and the maximum of each level is computed.
// Finally the addition of the maximum is computed to obtain the estimation
// of the parallel cost, and the final addition to obtain the sequential cost.
// These additions accumulate the counters by subsets (n1,n2) and (n3,n4),
// and the results are returned in the parameters auxS and auxP.
void PrintBlockTime (int typeP, ptr_ILU0Factor vFact, int nTim1, int nTim2,
											matDoubles mSumMax, int hmax, matDoubles *mSumThr, int nthr,
									int n1, int n2, int n3, int n4, double *auxS, double *auxP) {
	// Definition of the global vectors
	int *hgthtab, *ownrtab;
	// Definition of the local vectors and variables
	int i, j, k;
	double total_seq = 0.0, total_par = 0.0;

	// Parameters validation
	if ((typeP < 0) || (typeP > 4) || (vFact == NULL) || (vFact->mTab == NULL) || 
			(nTim1 < 0) || (nTim1 >= SIZE_TIM_LOC) || (nTim2 < 0) || (nTim2 > SIZE_TIM_LOC) || 
			(mSumMax == NULL) || (hmax < 0) || (mSumThr == NULL) || (nthr < 0) || 
			(n1 < 0) || (n1 >= SIZE_TIM_LOC) || (n2 < 0) || (n2 > SIZE_TIM_LOC) || 
			(auxS == NULL) || (auxP == NULL)) {
		printf ("Incorrect parameters in PrintBlockTime (%d,%d,%d,%d,%d,%d,%d)\n", 
										typeP, nTim1, nTim2, hmax, nthr, n1, n2); 
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	hgthtab = vFact->mTab[HGTH]; ownrtab = vFact->mTab[OWNR];
	// Initialize the auxiliar matrix
	InitDoubles (mSumMax[0], hmax*SIZE_TIM_LOC, 0.0, 0.0);
	InitDoubles (mSumMax[hmax],   SIZE_TIM_LOC, 0.0, 0.0);
	// Print the auxiliar vectors
	if (typeP < 1) PrintInts (ownrtab, vFact->dimL); 
	if (typeP < 1) PrintInts (hgthtab, vFact->dimL);
	// Print the information for each node and
	// accumulate the times for each thread
	for (j=0; j<vFact->dimL; j++) {
		if (typeP < 1) printf ("Node %d\n", j);
		for (i=nTim1; i<nTim2; i++) {
			if (typeP < 1) printf ("%20.10e ", vFact[j].tLoc[0][i]);
			if (ownrtab[j] < 0)
        mSumThr[0][hgthtab[j]-1][i] += vFact[j].tLoc[0][i];
      else
        mSumThr[ownrtab[j]][hgthtab[j]-1][i] += vFact[j].tLoc[0][i];
		}
		if (typeP < 1) printf ("\n");
	}
	// Print the auxiliar vectors
	if (typeP < 2) PrintInts (ownrtab, vFact->dimL); 
	if (typeP < 2) PrintInts (hgthtab, vFact->dimL);
	// Print the information for each Thread and
	// accumulate the times to know the maximum and the addition
	for (k=0; k<nthr; k++) {
		if (typeP < 2) printf ("Thread %d\n", k);
		for (j=0; j<hmax; j++) {
			for (i=nTim1; i<nTim2; i++) {
				if (typeP < 2) printf ("%20.10e ", mSumThr[k][j][i]);
				if (mSumMax[j][i-nTim1] < mSumThr[k][j][i])
					mSumMax[j][i-nTim1] = mSumThr[k][j][i];
				mSumMax[hmax][i-nTim1] += mSumThr[k][j][i];
			}
			if (typeP < 2) printf ("\n");
		}
	}
	// Compute the sequential cost as the addition of the previously added values
	total_seq = AddDoubles (mSumMax[hmax], n2-n1);
	if (n3 >= 0) total_seq += AddDoubles (mSumMax[hmax]+(n3-n1), n4-n3);
	if (typeP < 3) printf ("Total Seq Time = %20.10e = ", total_seq);
	if (typeP < 3) PrintDoubles (mSumMax[hmax], nTim2-nTim1);
	CopyDoubles (mSumMax[0], mSumMax[hmax], nTim2-nTim1);
	for (j=1; j<hmax; j++)
		for (i=0; i<nTim2-nTim1; i++) 
			mSumMax[hmax][i] += mSumMax[j][i];
	// Compute the sequential cost as the addition of the maximum of each level
	total_par = AddDoubles (mSumMax[hmax], n2-n1);
	if (n3 >= 0) total_par += AddDoubles (mSumMax[hmax]+(n3-n1), n4-n3);
	if (typeP < 3) printf ("Total Par Time = %20.10e = ", total_par);
	if (typeP < 3) PrintDoubles (mSumMax[hmax], nTim2-nTim1);
	if (typeP < 3) printf ("Increment (estimation) = %20.10e", total_seq / total_par);
	// Returns the sequential and parallel cost.
	*auxS = total_seq; *auxP = total_par;
}

// This routine prints the information of the counters, acording to the value of typeP:
// * < 1, it prints all the information, from the nodes to the final summarize.
// * < 2, it prints the information, from the threads to the final summarize.
// * < 3, it only prints the summarize.
void PrintTimesILU0Factor (int typeP, ptr_ILU0Factor vFact, int itrEnd) {
	// Definition of the global vectors
	int *hgthtab, *ownrtab;
	// Definition of the local vectors and variables
	int i, hmax, nthr;
	double total_seq, total_par, auxP = 0.0, auxS = 0.0;
	matDoubles mSumMax = NULL, *mSumThr = NULL;

	// Parameters validation
	if ((typeP < 0) || (typeP > 3) || (vFact == NULL) || (vFact->mTab == NULL)) {
		printf ("Incorrect parameters in PrintTimesILU0Factor (%d)\n", typeP); 
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	hgthtab = vFact->mTab[HGTH]; ownrtab = vFact->mTab[OWNR];
	// Compute the sizes of the rows of the auxiliar matrices
	hmax = MaxInts (hgthtab, vFact->dimL);
 	nthr = MaxInts (ownrtab, vFact->dimL) + 1;
	// Create the auxiliar matrix
	CreateMatrixDoubles (&mSumMax, hmax+1, SIZE_TIM_LOC);
	InitDoubles (mSumMax[0], (hmax+1)*SIZE_TIM_LOC, 0.0, 0.0);
	// Create the vector of auxiliar matrices
	mSumThr = (matDoubles *) malloc (nthr * sizeof(matDoubles));
	for (i=0; i<nthr; i++)  { 
		CreateMatrixDoubles (mSumThr+i, hmax, SIZE_TIM_LOC);
		InitDoubles (mSumThr[i][0], SIZE_TIM_LOC*hmax, 0.0, 0.0);
	}

	// Print the Global times
	printf ("Global Times = ");
	for (i=0; i<SIZE_TIM_GLB; i++)
		printf ("%20.10e ", vFact->tGlb[0][i]);
	printf ("\n\n");

	// Print the times related to the computation of the preconditioner
	printf ("Preconditioner Times\n");
	PrintBlockTime (typeP, vFact, TLBLDP, TLAMGP+1, mSumMax, hmax, mSumThr, nthr, 
							TLBLDP, TLFCTP+1, -1, -1, &total_seq, &total_par);
	if (typeP < 3) printf (" , ");
	printf ("Increment (real) = %20.10e , ", total_seq / vFact->tGlb[0][TGPPRC]);
	printf ("Percentage = %20.10e \n\n", total_par/vFact->tGlb[0][TGPPRC]);

	// Print the times related to the application of the preconditioner
	if (typeP < 3) printf ("Resolution/Transform Times (%d)\n", itrEnd);
	PrintBlockTime (typeP, vFact, TLBLDU, TLSOLD+1, mSumMax, hmax, mSumThr, nthr, 
								TLBLDU, TLRESU+1, TLBLDD, TLRESD+1, &auxS, &auxP);
	if (typeP < 3) printf ("\n\n");

	// Print the times related to the PCG
	printf ("PCG Times (%d,%d,%d)\n", nthr, vFact->dimL, itrEnd);
	PrintBlockTime (typeP, vFact, TLPROD, TLVEC5+1, mSumMax, hmax, mSumThr, nthr, 
								TLPROD, TLVEC5+1, -1, -1, &total_seq, &total_par);
	if (typeP < 3) printf ("\n");
	printf ("Increment (real) = %20.10e , ", (auxS+total_seq) / vFact->tGlb[0][TGPPCG]);
	printf ("Percentage = %20.10e\n\n", (auxP+total_par)/vFact->tGlb[0][TGPPCG]);

	// Remove the auxiliar matrix
	RemoveMatrixDoubles (&mSumMax);
	// Remove the vector of auxiliar matrices
	for (i=nthr-1; i>=0; i--) RemoveMatrixDoubles (mSumThr+i);
	free (mSumThr); mSumThr = NULL;
}

/****************************************************************************/

