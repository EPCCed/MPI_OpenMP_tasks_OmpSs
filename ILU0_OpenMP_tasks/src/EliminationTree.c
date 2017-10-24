#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "reloj.h"
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "EliminationTree.h"

// #define ADJUST_TESTVECTOR 1
// #define ADD_TESTVECTORS 1
// #define PRINT_BUILDCOPY_VECTORS 1

/*********************************************************************************/

// Compute, recursively, the height of each node of the tree, where
// * chldtab[i] includes the head of the list of the children of node i-th
// * brthtab includes the rest of the lists which are headed in chldtab
// * hgthtab[i] includes the heigth of the node i-th
void ComputeHeigthNodes (int *chldtab, int *brthtab, int *hgthtab, int root) {
	// Definition of the local variables
	int val, node;

	// Parameters validation
	if ((chldtab == NULL) || (brthtab == NULL) || (hgthtab == NULL) || (root < 0)) {
		printf ("Incorrect parameters in ComputeHeigthNodes (%d)\n", root);
		PrintTrace (); exit (-1);
	}
	// val --> the height of the children is equal to 1 plus the height of the root.
	// node --> begins with the first element of the list of the children.
	val = hgthtab[root] + 1; node = chldtab[root];
	while (node != -1) {
		hgthtab[node] = val;  // Fix the height of node
		ComputeHeigthNodes (chldtab, brthtab, hgthtab, node);
		node = brthtab[node]; // Move to the next element of the list
	}
}

// Compute the sizes which are required to create the structures of the node tsk
// * treetab[i] includes the father of the node i-th 
// * sizetab[i] includes the number of rows in the node i-th 
// * wgthtab[i] includes the number of nonzeros in the node i-th 
// * dimM is the number of rows to be computed
// * nnzM is the number of nonzeros to be computed
// * nlev is the number of nodes to arrive to the root
void ComputeSizesNodeEliminationTree (int tsk, int *treetab, int *sizetab, int *wgthtab,
																			int *dimM, int *nnzM, int *nlev) {
	// Definition of the local vectors and variables
	int dim = 0, nnz = 0, lev = 0, node = tsk;

	// Parameters validation
	if ((treetab == NULL) || (sizetab == NULL) || (wgthtab == NULL) || (tsk < 0)) {
		printf ("Incorrect parameters in ComputeSizesNodeEliminationTree (%d)\n", tsk);
		PrintTrace (); exit (-1);
	}
	// From the node tsk to the root of the tree, 
	// the number of rows and the nonzeros are acumulated
	do {
		dim += sizetab[node]; nnz += wgthtab[node];
		lev++; node = treetab[node];
	} while (node != -1);
	*dimM = dim; *nnzM = nnz; *nlev = lev;
}

// From the treetab, this routine computes other vectors related with the tree
// * chldtab[i] includes the head of the list in which the children of node i-th are
// * brthtab includes the rest of the lists which are headed in chldtab
// * nmchtab[i] includes the size of these lists
// * hgthtab[i] includes the height of each node
// The routine returns the number of nodes of the elimination tree
int ComputeEliminationTreeVectors (int *treetab, int *chldtab, int *nmchtab, int *brthtab, int *hgthtab, int size) {
	// Definition of the local vectors and variables
	int i, node = 0;

	// Parameters validation
	if ((treetab == NULL) || (chldtab == NULL) || (nmchtab == NULL) || 
			(brthtab == NULL) || (hgthtab == NULL) || (size < 1)) {
		printf ("Incorrect parameters in ComputeEliminationTreeVectors (%d)\n", size);
		PrintTrace (); exit (-1);
	}
	// Initialize the vectors chldtab and brthtab
	InitInts (chldtab, size, -1, 0); InitInts (brthtab, size, -1, 0);
	// Visit of the nodes of the vector treetab, until the root is reached
	i=0;
	while ((i < size) && (node >= 0)) {
		node = treetab[i]; 
		if (node >= 0) {  // the node is not the root
			if (chldtab[node] < 0) { // Create a new list of children related to node
				chldtab[node] = i; nmchtab[node] = 1;
			} else { // Insert i as the header of the list related to node
				brthtab[i] = chldtab[node];
				chldtab[node] = i; nmchtab[node]++;
			}	
		}
		i++;
	}
	// Compute the vector hgthtab
	hgthtab[i-1] = 1;
	ComputeHeigthNodes (chldtab, brthtab, hgthtab, i-1);

	// Return the number of nodes of the tree
	return i;
}

/*********************************************************************************/

// This routine computes the number of nonzeros of each row block of the
// sparse matrix spr. The block structure in defined by the permutation
// perm and the auxiliar vectors rangtab and sizetab, while the result
// is incoporated to the auxiliar vector wgthtab. The size of all the
// auxiliar vectors is dimL.
// The parameter index indicates if 0-indexing or 1-indexing is used.
void ComputeWeightNodes (SparseMatrix spr, int index, int *perm, int *rangtab, int *sizetab, 
													int *wgthtab, int dimL) {
	// Definition of the local vectors and variables
	int i;

	// Parameters validation
	if ((spr.dim1 < 1) || (spr.vptr == NULL) || (spr.vpos == NULL) || (spr.vval == NULL) || 
			(rangtab == NULL) || (sizetab == NULL) || (wgthtab == NULL) || (dimL < 1)) {
		printf ("Incorrect parameters in ComputeEliminationTreeVectors (%d)\n", dimL);
		PrintTrace (); exit (-1);
	}
	// Compute the maximum number of nonzeros of each node (wgthtab vector)
	TransformHeadertoLength (spr.vptr, spr.dim1);
	for (i=0; i<dimL; i++) {
		wgthtab[i] = AddPermuteInts (spr.vptr+1, perm+rangtab[i], index, sizetab[i]);
	}
	TransformLengthtoHeader (spr.vptr, spr.dim1);
}

// Scale the differents blocks of the sparse matrix to transform to unconsistent:
// * task define the node in the elimination tree
// * treetab and nmchtab are the required vector related to the elimination tree
// * head is the vector which marks the begin and the end of each level (CSR way)
// * dimN and nlev are the dimensions of the previous vectors
// The parameter index indicates if 0-indexing or 1-indexing is used.
void AdjustLeafEliminationTree (SparseMatrix spr, int index, int task, int *treetab, 
																int *nmchtab, int dimN, int *head, double *divs, int nlev) { 
	// Definition of the local vectors and variables
	int i, j, k, lev; 
	int *pi2 = NULL, *pi3 = NULL, *pi4 = NULL, *pi5 = NULL, *pi6 = NULL;
	double *pd6 = NULL;

	// Parameters validation
	if ((spr.dim1 < 1) || (spr.vptr == NULL) || (spr.vpos == NULL) || (spr.vval == NULL) || 
			((index & (~1)) != 0) || (treetab == NULL) || (nmchtab == NULL) || (head == NULL) || 
			(divs == NULL) || (task < 0) || (dimN < 1) || (nlev < 1)) {
		printf ("Incorrect parameters in AdjustLeafEliminationTree (%d,%d,%d,%d)\n", 
							index, task, dimN, nlev);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	pi2 = head; pi3 = pi2+1;
	// Compute the divisor related to each level
	divs[0] = 1.0; lev = treetab[task];
	for (i=1; i<nlev; i++) {
		divs[i] = divs[i-1] / nmchtab[lev];
		lev = treetab[lev];
	}
	// Scale the blocks using the divisors
	pi2 = (pi3++); 
	for (i=1; i<nlev; i++) { // Avoid the first block
		pi4 = spr.vptr + (*pi2); pi5 = pi4+1;
		for (j=(*pi2); j<(*pi3); j++) { // head[i]-head[i+1]
			pi6 = spr.vpos - index + (*pi4); pd6 = spr.vval - index + (*pi4);
			for (k=(*pi4); k<(*pi5); k++) { // vptr[i]-vptr[i+1]
				// Locate the divisor to be applied
				lev = 0; 
				while ((lev < i) && ((*pi6 - index) >= head[lev+1])) lev++; 
				// if the divisor is bigger than 1.0
				if (lev > 0) (*pd6) *= divs[lev];
				pi6++; pd6++;
			}
			pi4 = (pi5++);
		}
		pi2 = (pi3++);
	}
	// Remove the negative indices (unnecessary in this case)
	PermuteColsWithNegSparseMatrix (spr, index, NULL);
}

// Remove the contraint blocks for the leaves which are in the nonminimum block.
// * task define the node in the elimination tree
// * treetab, chldtab and brthtab are the required vector related to the elimination tree
// * head is the vector which marks the begin and the end of each level (CSR way)
// * dimN and nlev are the dimensions of the previous vectors
// The parameter index indicates if 0-indexing or 1-indexing is used.
void AdjustLeafPCGEliminationTree (SparseMatrix spr, int index, int task, int *treetab,
																		int *chldtab, int *brthtab, int dimN, int *head, int nlev) {
	// Definition of the local vectors and variables
	int i, j, k, lev, node;
	int *pi2 = NULL, *pi3 = NULL, *pi4 = NULL, *pi5 = NULL, *pi6 = NULL, *minmtab = NULL;

	// Parameters validation
	if ((spr.dim1 < 1) || (spr.vptr == NULL) || (spr.vpos == NULL) || (spr.vval == NULL) ||
			((index & (~1)) != 0) || (treetab == NULL) || (chldtab == NULL) || (brthtab == NULL) ||
			(head == NULL) || (task < 0) || (dimN < 1) || (nlev < 1)) {
		printf ("Incorrect parameters in AdjustLeafEliminationTree (%d,%d,%d,%d)\n",
							index, task, dimN, nlev);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	CreateInts (&minmtab, nlev);
	// Compute the minimum task to each level
	minmtab[0] = task; node = task;
	for (i=1; i<nlev; i++) {
		j = chldtab[treetab[node]]; minmtab[i] = j;
		while (brthtab[j] != -1) {
			if (minmtab[i] > brthtab[j])
				minmtab[i] = brthtab[j];
			j = brthtab[j];
		}
		node = treetab[node];
	}
	// Remove the nonzeros in the constraint area for the nonminimum task
	pi2 = (head+1); pi3 = pi2+1;
	for (i=1; i<nlev; i++) { // Avoid the first block
		pi4 = spr.vptr + (*pi2); pi5 = pi4+1;
		for (j=(*pi2); j<(*pi3); j++) { // head[i]-head[i+1]
			pi6 = spr.vpos - index + (*pi4);
			for (k=(*pi4); k<(*pi5); k++) { // vptr[i]-vptr[i+1]
				// Locate the divisor to be applied
				lev = 0; node = task;
				while ((lev < i) && (node == minmtab[lev]) && ((*pi6 - index) >= head[lev+1])){
					if (lev > 0) node = treetab[node];
					lev++;
				}
				// if the divisor is bigger than 1.0
				if ((lev <= i) && (node != minmtab[lev])) (*pi6) = -1;
				pi6++;
			}
			pi4 = (pi5++);
		}
		pi2 = (pi3++);
	}
	// Remove the auxiliar vector
	RemoveInts (&minmtab);
	// Remove the negative indices
	PermuteColsWithNegSparseMatrix (spr, index, NULL);
}

// From the original data (spr,perm) and the elimination tree information (treetab,rangtab),
// the routine builds the data related to a leaf:
// * sprM is the sparse matrix, and permM is the corresponding permutation.
// * headL is the vector which marks the begin and the end of each level (CSR way)
// The parameters index and indexM indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices spr and sprM.
void BuildLeafFromMatrix (SparseMatrix spr, int index, int *perm, int *treetab, int *rangtab, 
													int leaf, SparseMatrix sprM, int indexM, int *permM, int *headL) {
	// Definition of the global and local vectors and variables
	int i, nlev = 0, node, dim = 0, row, sizerow, indexIO = indexM - index;
	int *ph, *pp, *pp1, *pi1, *pp2, *pi2;
	double *pd1, *pd2;

	// Parameters validation
	if ((spr.dim1 < 1) || (spr.vptr == NULL) || (spr.vpos == NULL) || (spr.vval == NULL) || 
			((index & (~1)) != 0) || (sprM.vptr == NULL) || (sprM.vpos == NULL) || (sprM.vval == NULL) || 
			((indexM & (~1)) != 0) || (perm == NULL) || (treetab == NULL) || (rangtab == NULL) || 
			(permM == NULL) || (headL == NULL) || (leaf < 0)) {
		printf ("Incorrect parameters in BuildLeafFromMatrix (%d,%d,%d)\n", index, leaf, indexM);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	ph = headL+1; pp = permM;
	pp1 = spr.vptr-index; pi1 = spr.vpos-index; pd1 = spr.vval-index; 
	pp2 = sprM.vptr+1; pi2 = sprM.vpos; pd2 = sprM.vval;
	// The elimination tree is traverse from the leaf to the root, 
	node = leaf;
	do {
		for (i=rangtab[node]; i<rangtab[node+1]; i++) { // All the rows of the node
			row = perm[i]; *(pp++) = row + indexIO;
			*(pp2++) = sizerow = pp1[row+1] - pp1[row];
			CopyShiftInts (pi1 + pp1[row], pi2, sizerow, indexIO); pi2 += sizerow;
			CopyDoubles (pd1 + pp1[row], pd2, sizerow); pd2 += sizerow;
		}
		*ph = rangtab[node+1]-rangtab[node]; dim += *(ph++); nlev++;
		node = treetab[node];
	} while (node != -1);
	// Transform the lenght vectors to the CSR way
	*(sprM.vptr) = indexM; TransformLengthtoHeader (sprM.vptr, dim);
	*(headL) = 0; TransformLengthtoHeader (headL, nlev);
}

// From the matrix spr, the routine computes the leaf related to the node task,
// and the resulted matrix is included in the node task
// The parameter indexM indicates if 0-indexing or 1-indexing is used to extract the factor
// If the parameter tvG is not null, it contains the test vector to be used in 
// the preconditioner computation.
void BuildLeafFromMatrixF (SparseMatrix spr, int indexM, double *tvG, ptr_ILU0Factor vFact, int task) {
	// Definition of the global vectors and variables
	int dim, dimL, index;
	int *perm, *treetab, *rangtab, *sizetab, *nmchtab, *wgthtab;
	// Definition of the local vectors and variables
	int i, dimM, dimD, nnzM, nlev;
	int *permM = NULL, *ipermM = NULL, *headL = NULL;
	double *diag = NULL, *diag1 = NULL, *diag2 = NULL, *divs = NULL, *tv = NULL;
	matDoubles mDia;
	SparseMatrix sprM;

	// Parameters validation
	if ((spr.dim1 < 1) || (spr.vptr == NULL) || (spr.vpos == NULL) || (spr.vval == NULL) || 
			((indexM & (~1)) != 0) || (vFact == NULL) || ((vFact->indP & (~1)) != 0) || 
			(vFact->perm == NULL) || (vFact->mTab == NULL) || (vFact->dimL < 1) || (task < 0)) {
		printf ("Incorrect parameters in BuildLeafFromMatrixF (%d,%d)\n", indexM, task);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	dim = spr.dim1; dimL = vFact[task].dimL; index = vFact->indP;
	perm = vFact[task].perm; treetab = vFact[task].mTab[TREE]; rangtab = vFact[task].mTab[RANG];
	sizetab = vFact[task].mTab[SIZE]; nmchtab = vFact[task].mTab[NMCH]; wgthtab = vFact[task].mTab[WGTH];
	// Calculate the properties of the node
	ComputeSizesNodeEliminationTree (task, treetab, sizetab, wgthtab, &dimM, &nnzM, &nlev);
	// Create the local structures
	CreateSparseMatrix (&sprM, indexM, dimM, dimM, nnzM, 0);
	CreateInts (&permM, dimM+nlev+1); headL = permM + dimM;
	dimD = dimM - sizetab[task]; 
	CreateMatrixDoubles (&mDia, SIZE_MAT_DIA, dimD);
	diag = mDia[DIAG]; diag1 = mDia[DIA1]; diag2 = mDia[DIA2];
	CreateDoubles (&divs , nlev); CreateInts (&ipermM, dim); 
	// Fill the local structures from the shared data
	BuildLeafFromMatrix (spr, index, perm, treetab, rangtab, task, sprM, indexM, permM, headL);
	// Eliminate incorrect elements of Sparse Matrices
	for (i=0; i<dim; i++) ipermM[i] = -1;
	ComputeInvPermutation (permM, ipermM, indexM, dimM);
	PermuteColsWithNegSparseMatrix (sprM, indexM, ipermM);
	// Readjust the size of the local data
	ReallocSparseMatrix (&sprM);
	// Scale the replicate blocks, getting the diagonal before and after the operation
	GetDiagonalSparseMatrixDspls (sprM, indexM, diag1, sizetab[task]);
	AdjustLeafEliminationTree (sprM, indexM, task, treetab, nmchtab, dimL, headL, divs, nlev);
	GetDiagonalSparseMatrixDspls (sprM, indexM, diag2, sizetab[task]);
	// Compute the vector diag as (diag1-diag2)
	CopyDoubles (diag1, diag, dimD); AxpyDoubles (-1.0, diag2, diag, dimD);
	// Free unuseful data
	RemoveInts (&ipermM);
	// Computation of the test vector
	if (tvG != NULL) {
		CreateDoubles (&tv, dimM); CopyPermuteDoubles (tvG, permM, indexM, tv, dimM);
#ifdef ADJUST_TESTVECTOR
		AdjustVector (tv, divs, headL, nlev);
#endif
	}

	// Repair the structure
	vFact[task].sprM = sprM  ; vFact[task].dimM  = dimM ; vFact[task].permM = permM; 
	vFact[task].indM = indexM; 												   vFact[task].mDia  = mDia  ; 
	vFact[task].nlev = nlev  ; vFact[task].headL = headL; vFact[task].divs  = divs ;
}
#define NEW_ADJUST_PCG 1  //SE COMENTA SI QUEREMOS MISMA MATRIZ QUE PRECOND
// From the matrix spr, the routine computes the PCG leaf related to the node task,
// and the resulted matrix is included in the node task
// The parameter indexM indicates if 0-indexing or 1-indexing is used to extract the factor
void BuildLeafFromMatrixPCG (SparseMatrix spr, int indexR, ptr_ILU0Factor vFact, int task) {
	// Definition of the global vectors and variables
	int dim, dimL, index;
#ifdef NEW_ADJUST_PCG
	int *perm, *treetab, *rangtab, *sizetab, *wgthtab, *brthtab, *chldtab;
#else
	int *perm, *treetab, *rangtab, *sizetab, *wgthtab;
#endif
	// Definition of the local vectors and variables
	int i, dimM, nnzR, nlev;
	int *permR = NULL, *ipermR = NULL, *headL = NULL;
	SparseMatrix sprR;

	// Parameters validation
	if ((spr.dim1 < 1) || (spr.vptr == NULL) || (spr.vpos == NULL) || (spr.vval == NULL) || 
			((indexR & (~1)) != 0) || (vFact == NULL) || ((vFact->indP & (~1)) != 0) || 
			(vFact->perm == NULL) || (vFact->mTab == NULL) || (vFact->dimL < 1) || (task < 0)) {
		printf ("Incorrect parameters in BuildLeafFromMatrixF (%d,%d)\n", indexR, task);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	dim = spr.dim1; dimL = vFact[task].dimL; index = vFact->indP;
	perm = vFact[task].perm; treetab = vFact[task].mTab[TREE]; rangtab = vFact[task].mTab[RANG];
	sizetab = vFact[task].mTab[SIZE]; wgthtab = vFact[task].mTab[WGTH];
#ifdef NEW_ADJUST_PCG
	chldtab = vFact[task].mTab[CHLD]; brthtab = vFact[task].mTab[BRTH];
#endif
	// Calculate the properties of the node
	ComputeSizesNodeEliminationTree (task, treetab, sizetab, wgthtab, &dimM, &nnzR, &nlev);
	// Create the local structures
	CreateSparseMatrix (&sprR, indexR, dimM, dimM, nnzR, 0);
	permR = vFact[task].permM; headL = vFact[task].headL;
	CreateInts (&ipermR, dim); 
	// Fill the local structures from the shared data
	BuildLeafFromMatrix (spr, index, perm, treetab, rangtab, task, sprR, indexR, permR, headL);
	// Eliminate incorrect elements of Sparse Matrices
	for (i=0; i<dim; i++) ipermR[i] = -1;
	ComputeInvPermutation (permR, ipermR, indexR, dimM);
	PermuteColsWithNegSparseMatrix (sprR, indexR, ipermR);
#ifdef NEW_ADJUST_PCG
	// Remove the replicate blocks, such that only a copy exists
	AdjustLeafPCGEliminationTree (sprR, indexR, task, treetab, chldtab, brthtab, dimL, headL, nlev);
#else
	// Scale the replicate blocks, getting the diagonal before and after the operation
	AdjustLeafEliminationTree (sprR, indexR, task, treetab, nmchtab, dimL, headL, divs, nlev);
#endif
	// Readjust the size of the local data
	ReallocSparseMatrix (&sprR);
	// Free unuseful data
	RemoveInts (&ipermR); 

	// Repair the structure
	vFact[task].sprR = sprR; vFact[task].indR = indexR;
}

//#define JOSE_FLOPS 1
#define EASY_SORT 1

// This routine computes the matrix to be processed, adding the matrices computed 
// by the children of the task and the result is included in the position task of vFact. 
// The parameter indexM indicates if 0-indexing or 1-indexing is used to accumulate the matrices.
// This routine returns the cost of the addition of the matrices of the children
double BuildNodeFromChildren (int indexM, ptr_ILU0Factor vFact, int task) {
	// Definition of the global vectors and variables
	int *treetab, *sizetab, *chldtab, *brthtab, *wgthtab;
	// Definition of the local vectors and variables
	int i, dimM, dimD, nnzM, nlev;
	int indexF1 = 0, indexF2 = 0;
	int node, node1, node2, dim1, dim2R, nnz2R, dim2L, nnz2L, dim3, nmch = 0;
	int *perm1 = NULL, *iperm1 = NULL, *perm2 = NULL, *iperm2 = NULL, *headL = NULL;
	SparseMatrix sprM, sprM1, sprM2;
	int *permM = NULL, *ipermM = NULL;
	double *diag = NULL, *diag1 = NULL, *diag2 = NULL, *divs = NULL;
	double *tv1 = NULL, *tv2 = NULL, *tv = NULL;
	matDoubles mDia;
	double te = 0.0, tu = 0.0, te1, te2, tu1, tu2;
#ifdef EASY_SORT
	int rangs1[4], shfts1[4], rangs2[4], shfts2[4]; 
#endif
#ifdef JOSE_FLOPS
	long int Jcont = 0;
#endif
	
	// Parameters validation
	if (((indexM & (~1)) != 0) || (vFact == NULL) || (vFact->mTab == NULL) || (task < 0)) {
		printf ("Incorrect parameters in BuildNodeFromChildren (%d,%d,%d,%d)\n", 
							indexM, task, (vFact == NULL), (vFact->mTab == NULL));
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	treetab = vFact[task].mTab[TREE];
	sizetab = vFact[task].mTab[SIZE]; chldtab = vFact[task].mTab[CHLD];
	brthtab = vFact[task].mTab[BRTH]; wgthtab = vFact[task].mTab[WGTH];
	// Calculate the properties of the node
	ComputeSizesNodeEliminationTree (task, treetab, sizetab, wgthtab, &dimM, &nnzM, &nlev);
	// Create new structures
	dim1 = sizetab[task]; dim3 = dimM - dim1; dimD = dim3;
	CreateMatrixDoubles (&mDia, SIZE_MAT_DIA, dimD);
	diag = mDia[DIAG]; diag1 = mDia[DIA1]; diag2 = mDia[DIA2];
	CreateDoubles (&divs, nlev); 

	// Acumulate the matrices computed by the children
	node = task; node1 = chldtab[node]; node2 = brthtab[node1]; 
	while (node2 != -1) {

		// Select the operators in the acumulation
		sprM1   = (nmch == 0)? vFact[node1].sprF: sprM  ; sprM2   = vFact[node2].sprF; 
		indexF1 = (nmch == 0)? vFact[node1].indF: indexM; indexF2 = vFact[node2].indF; 
		// Data validation
		if ((sprM1.dim1 < 1) || (sprM1.vptr == NULL) || (sprM1.vpos == NULL) || (sprM1.vval == NULL) || 
				(sprM2.dim1 < 1) || (sprM2.vptr == NULL) || (sprM2.vpos == NULL) || (sprM2.vval == NULL) ||
				((indexF1 & (~1)) != 0) || ((indexF2 & (~1)) != 0)) {
			printf ("Incorrect parameters in BuildNodeFromChildren (%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)\n", 
								task, node1, node2, indexF1, indexF2, 
								(sprM1.dim1 < 1), (sprM1.vptr == NULL), (sprM1.vpos == NULL), (sprM1.vval == NULL),
        				(sprM2.dim1 < 1), (sprM2.vptr == NULL), (sprM2.vpos == NULL), (sprM2.vval == NULL));
			PrintTrace (); exit (-1);
		}
		// Fill the diag's and perm's vectors
		if (nmch == 0) {
			if (dim3 > 0) {
				CopyDoubles (vFact[node1].mDia[DIA1]+dim1, diag1, dimD);
				CopyDoubles (vFact[node1].mDia[DIA2]+dim1, diag2, dimD);
			}
			permM = ipermM = NULL;
		}
		// Compute the size of the new matrix and create this matrix
		dim2R = sprM1.dim1 - dimM; nnz2R = sprM1.vptr[sprM1.dim1]-indexF1;
		dim2L = sprM2.dim1 - dimM; nnz2L = sprM2.vptr[sprM2.dim1]-indexF2;
		CreateSparseMatrix (&sprM, indexM, dimM+dim2R+dim2L, dimM+dim2R+dim2L, nnz2R+nnz2L, 0); 
		// Creation of the test vector
		if (tv1 != NULL) CreateDoubles (&tv, dimM+dim2R+dim2L); 
		// Define the permutation of the right matrix or the previous addition
		// and the corresponding test vector
		CreateInts (&perm1, sprM1.dim1); CreateInts (&iperm1, sprM1.dim1+dim2L); 
		for (i=0; i<dim2R; i++) perm1[i] = i+indexF1;
		for (i=dim2R; i<dim2R+dim1; i++) perm1[i] = i+dim2L+indexF1;
		if (tv1 != NULL) { 
			CopyDoubles (tv1, tv, dim2R); CopyDoubles (tv1+dim2R, tv+dim2R+dim2L, dim1);
		}
		for (i=dim1+dim2R; i<dim3+dim2R+dim1; i++) perm1[i] = i+dim2L+indexF1;
		if (tv1 != NULL) CopyDoubles (tv1+dim2R+dim1, tv+dim2R+dim2L+dim1, dim3);
		InitInts (iperm1, sprM1.dim1+dim2L, -1, 0);
		ComputeInvPermutation (perm1, iperm1, indexF1, sprM1.dim1);
		// Define the permutation of the left matrix and the corresponding test vector
		CreateInts (&perm2, sprM2.dim1); CreateInts (&iperm2, sprM2.dim1+dim2R); 
		for (i=0; i<dim2L+dim1; i++) perm2[i] = i+dim2R+indexF2;
#ifdef ADD_TESTVECTORS
		if (tv2 != NULL) {
	#ifdef JOSE_FLOPS
			Jcont += dim1;
	#endif
			CopyDoubles (tv2, tv+dim2R, dim2L); AxpyDoubles (1.0, tv2, tv+dim2R+dim2L, dim1);
		}
#else
		if (tv2 != NULL) CopyDoubles (tv2, tv+dim2R, dim2L+dim1);
#endif
		for (i=dim1+dim2L; i<dim3+dim2L+dim1; i++) perm2[i] = i+dim2R+indexF2;
#ifdef ADD_TESTVECTORS
	#ifdef JOSE_FLOPS
		if (tv2 != NULL) Jcont += dim3;
	#endif
		if (tv2 != NULL) AxpyDoubles (1.0, tv2+dim2L+dim1, tv+dim2R+dim2L+dim1, dim3);
#else
		if (tv2 != NULL) CopyDoubles (tv2+dim2L+dim1, tv+dim2R+dim2L+dim1, dim3);
#endif
		InitInts (iperm2, sprM2.dim1+dim2R, -1, 0);
		ComputeInvPermutation (perm2, iperm2, indexF2, sprM2.dim1);
		// Compute the addition
#ifdef JOSE_FLOPS
//		Jcont += IntersectSparseMatrices (sprM1, indexF1, perm1, iperm1, sprM2, 
//																			indexF2, perm2, iperm2, &sprM, indexM);
#endif
		reloj (&te1, &tu1);
		// Shift the vpos vectors
#ifdef EASY_SORT
		rangs1[0] = 0+indexF1; rangs1[1] = dim2R+indexF1; rangs1[2] = dim2R+dim1+indexF1;
		shfts1[0] = 0; shfts1[1] = dim2L; shfts1[2] = dim2L;
		CopyShiftRangsInts (sprM1.vpos, sprM1.vpos, nnz2R, rangs1, shfts1, 3);
		rangs2[0] = 0+indexF2; rangs2[1] = dim2L+indexF2; rangs2[2] = dim2L+dim1+indexF2;
		shfts2[0] = dim2R; shfts2[1] = dim2R; shfts2[2] = dim2R;
		CopyShiftRangsInts (sprM2.vpos, sprM2.vpos, nnz2L, rangs2, shfts2, 3);
#endif
		// Compute the addition
#ifdef EASY_SORT
		AddSparseMatrices (sprM1, indexF1, perm1, iperm1, sprM2, indexF2, perm2, iperm2, &sprM, indexM, 0);
#else
		AddSparseMatrices (sprM1, indexF1, perm1, iperm1, sprM2, indexF2, perm2, iperm2, &sprM, indexM, 1);
#endif
		// Unshift the vpos vectors
#ifdef EASY_SORT
		for (i=0; i<3; i++) { shfts1[i] = - shfts1[i]; shfts2[i] = - shfts2[i]; }
		CopyShiftRangsInts (sprM1.vpos, sprM1.vpos, nnz2R, rangs1, shfts1, 3);
		CopyShiftRangsInts (sprM2.vpos, sprM2.vpos, nnz2L, rangs2, shfts2, 3);
#endif
		reloj (&te2, &tu2); te += (te2 - te1); tu += (tu2 - tu1);
		ReallocSparseMatrix (&sprM);
		// Free the unseful data
		RemoveInts (&iperm2); RemoveInts (&perm2); RemoveInts (&iperm1); RemoveInts (&perm1);
		if (nmch == 0) {
			nmch = 2;
		} else {
			RemoveSparseMatrix (&sprM1); nmch++;
		}
		node1 = node2; node2 = brthtab[node1];
	}
#ifdef JOSE_FLOPS
	Jcont += dim3;
#endif
	// Scale the diag2
	ScaleDoubles (diag2, nmch*1.0, dim3);
	// Compute the vector diag as (diag1-diag2)
	CopyDoubles (diag1, diag, dimD); AxpyDoubles (-1.0, diag2, diag, dimD);
	// Computation of nlev, headL, divs
	CreateInts (&permM, nlev+1); headL = permM;
	for (i=2; i<=nlev; i++)
		headL[i] = vFact[node1].headL[i+1]-vFact[node1].headL[i];
	headL[0] = 0; headL[1] = sprM.dim1 - dim3;
	TransformLengthtoHeader (headL, nlev);
	for (i=0; i<nlev; i++) divs[i] = vFact[node1].divs[i+1] * nmch;

#ifdef JOSE_FLOPS
	printf ("Cflops_1 = %22ld\n", Jcont);
#endif
	// Repair the structure
	vFact[task].sprM = sprM  ; vFact[task].dimM  = dimM ; vFact[task].permM = permM; 
	vFact[task].indM = indexM; 												   vFact[task].mDia  = mDia  ;
	vFact[task].nlev = nlev  ; vFact[task].headL = headL; vFact[task].divs  = divs ;
	// Return the cost of the additions
	return te;
}

/*********************************************************************************/

// Scale the differents blocks of the vector to transform to unconsistent:
// * vecL is the vector to be adjusted
// * divs includes the divisors related to each level
// * head is the vector which marks the begin and the end of each level (CSR way)
// * nlev is the dimension of the previous vectors
void AdjustVector (double *vecL, double *divs, int *head, int nlev) { 
	// Definition of the local vectors and variables
	int i, *pi2, *pi3;

	// Parameters validation
	if ((vecL == NULL) || (divs == NULL) || (head == NULL) || (nlev < 1)) {
		printf ("Incorrect parameters in AdjustVector (%d)\n", nlev);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	pi2 = head; pi3 = pi2+1;
	// Scale the blocks using the divisors
	pi2 = (pi3++); 
	for (i=1; i<nlev; i++) { // Avoid the first block
		ScaleDoubles (vecL+(*pi2), divs[i], (*pi3)-(*pi2));
		pi2 = (pi3++);
	}
}

// From the original vector (vec) and the local permutation (permM), 
// the routine builds the local vector (vecL)
// The parameter index indicates if 0-indexing or 1-indexing is used.
void BuildLeafFromVector (double *vec, int *permM, int index, double *vecL, int dimM) {
	// Parameters validation
	if ((vec == NULL) || (permM == NULL) || (vecL == NULL) || (dimM < 1) ||
			((index & (~1)) != 0)) {
		printf ("Incorrect parameters in BuildLeafFromVector (%d,%d)\n", index, dimM);
		PrintTrace (); exit (-1);
	}
	CopyPermuteDoubles (vec, permM, index, vecL, dimM);
}

// From the original vector (vec), the routine extracts the components 
// corresponding to the permutation related to the node task, 
// and it stores them in the vector vid related to the node task.
// The vector will be scaled if the parameter adjustV is nonzero.
void BuildLeafFromVectorF (double *vec, ptr_ILU0Factor vFact, int task, 
														int vid, int adjustV) {
	// Definition of the global vectors and variables
	int dimM, indexM, nlev, *permM, *headL;
	double *divs, *ptr; 

	// Parameters validation
	if ((vec == NULL) || (vFact == NULL) || (task < 0) || (vid < 0) || 
			(vFact[task].permM == NULL) || (vFact[task].headL == NULL) || 
			(vFact[task].divs == NULL)) {
		printf ("Incorrect parameters in BuildLeafFromVectorF (%d,%d)\n", task, vid);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	dimM  = vFact[task].dimM ; indexM = vFact[task].indR ; nlev = vFact[task].nlev; 
	permM = vFact[task].permM; headL  = vFact[task].headL; divs = vFact[task].divs;
	// Create the local structures
	if (vFact[task].mVcL == NULL) {
		CreateMatrixDoubles (&(vFact[task].mVcL), SIZE_MAT_VCL, dimM);
		InitDoubles (*vFact[task].mVcL, SIZE_MAT_VCL*dimM, -2.0, 0.0);
		vFact[task].dimV = dimM; vFact[task].tskV = task;
	}
	// Determinate the vector on which the data will be stored
	ptr = IdentifyVectorResolution (vFact+task, vid);
	// Fill the local structures from the shared data
	BuildLeafFromVector (vec, permM, indexM, ptr, dimM);
	// Scale the vector
	if (adjustV) AdjustVector (ptr, divs, headL, nlev);
}

// This routine computes the vector to be processed, adding the vectors vidI computed 
// by the children of the task and the result is included in the vector vidO of vFact. 
// In the node task, the vector vidI is used as an auxiliar vector.
void BuildNodeFromChildrenV (ptr_ILU0Factor vFact, int task, int vidI, int vidO) {
	// Definition of the global vectors and variables
	int dimM, *sizetab, *chldtab, *nmchtab, *brthtab;
	// Definition of the local vectors and variables
	int node, node1, node2, nmch = 0, even_odd;
	int dimXX, dim1, dim2R, dim2L, dim3, dimF1, dimF2, dimM1, dimM2;
	double *pL, *pL1, *pL2, *pL3, *pR1, *pR2, *pR3, *ptrI, *ptrO;

	// Parameters validation
	if ((vFact == NULL) || (task < 0) || (vidI < 0) || (vidO < 0) || 
			(vFact->mTab == NULL)) {
		printf ("Incorrect parameters in BuildNodeFromChildrenV (%d,%d,%d)\n", 
								task, vidI, vidO);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	sizetab = vFact[task].mTab[SIZE]; chldtab = vFact[task].mTab[CHLD];
	nmchtab = vFact[task].mTab[NMCH]; brthtab = vFact[task].mTab[BRTH]; 
	dimM = vFact[task].dimM; dimXX =  vFact[task].dimX; 
	dim1 = sizetab[task]; dim3 = dimM - dim1;
	// Create new structures
	if (vFact[task].mVcL == NULL) {
		CreateMatrixDoubles (&(vFact[task].mVcL), SIZE_MAT_VCL, dimXX);
		InitDoubles (*vFact[task].mVcL, SIZE_MAT_VCL*dimXX, -2.0, 0.0);
		vFact[task].dimV = dimXX; vFact[task].tskV = task;
	}
	// Determinate the vectors vidI and vidO of node task 
	ptrI = IdentifyVectorResolution (vFact+task, vidI);
	ptrO = IdentifyVectorResolution (vFact+task, vidO);
	// Acumulate the result matrices of the children
	node = task; node1 = chldtab[node]; node2 = brthtab[node1];
	even_odd = (nmchtab[node] % 2);
	while (node2 != -1) {
		// Define the sizes and the sections of the vectors to be acumulated
		// Definition of the data for the first operand
		if (nmch == 0) { // If no accumulation has been made
			dimM1 = vFact[node1].dimX; dimF1 = vFact[node1].dimF;
			dim2R = dimF1 - dimM; 
			pR2 = IdentifyVectorResolution (vFact+node1, vidI)+(dimM1-dimF1); 
			pR1 = pR2 + dim2R; pR3 = pR1 + dim1;
		} else { // If some accumulation has been made
			dimF1 = dimXX; dim2R = dimF1 - dimM; 
			pR2 = (even_odd)? ptrO: ptrI; 
			pR1 = pR2 + dim2R; pR3 = pR1 + dim1 ;
		}
		// Definition of the data for the second operand
		dimM2 = vFact[node2].dimX; dimF2 = vFact[node2].dimF;
		dim2L = dimF2 - dimM; 
		pL2 = IdentifyVectorResolution (vFact+node2, vidI)+(dimM2-dimF2); 
		pL1 = pL2 + dim2L; pL3 = pL1 + dim1;
		// Select where acumulate the result and its size
		dimXX = dimM+dim2R+dim2L; 
		pL = (even_odd)? ptrI: ptrO; 
		// Acumulate the sectors of the vectors
		CopyDoubles (pR2, pL, dim2R); pL += dim2R;
		CopyDoubles (pL2, pL, dim2L); pL += dim2L;
		CopyDoubles (pR1, pL, dim1 ); AxpyDoubles (1.0, pL1, pL, dim1); pL += dim1;
		CopyDoubles (pR3, pL, dim3 ); AxpyDoubles (1.0, pL3, pL, dim3); pL += dim3;
		// Prepare the next iteration
		nmch = ((nmch == 0)? 1: nmch) + 1;
		node1 = node2; node2 = brthtab[node1]; even_odd = (1 - even_odd);
	}
}

// This routine computes the vector to be processed, copying the vector vidI computed
// by the father of node task on the vector vidO of node task.
void BuildNodeFromFatherV (ptr_ILU0Factor vFact, int task, int vidI, int vidO) {
	// Definition of the global vectors and variables
	int *treetab,  *sizetab, *chldtab, *brthtab;
	// Definition of the local vectors and variables
	int fthr, node; 
	int dim1, dim2R, dim3, dimXX;
	int dimF; 
	double *pL1, *pL2, *pL3, *pR1, *pR2, *pR3;
	int dimMF;
	double *ptrO;

	// Parameters validation
	if ((vFact == NULL) || (task < 0) || (vidI < 0) || (vidO < 0) || 
			(vFact->mTab == NULL)) {
		printf ("Incorrect parameters in BuildNodeFromFatherV (%d,%d,%d)\n", 
								task, vidI, vidO);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	treetab = vFact[task].mTab[TREE]; sizetab = vFact[task].mTab[SIZE];
	chldtab = vFact[task].mTab[CHLD]; brthtab = vFact[task].mTab[BRTH];
	// Calculate the properties of the node
	fthr = treetab[task]; 
	dimMF = vFact[fthr].dimM; dimXX = vFact[fthr].dimX;
	dim1 = sizetab[fthr]; dim3  = vFact[fthr].dimM - dim1; 
	// Define the sectors of the vector in the node fthr
	pL2 = IdentifyVectorResolution (vFact+fthr, vidI);
	pL3 = pL2 + (dimXX - dim3); pL1 = pL3 - dim1; 
	// Define where the sectors of bad pivots related to node task is
	// and its size
	node = chldtab[fthr]; 
	while (node != task) {
		dimF = vFact[node].dimF; dim2R = dimF - dimMF; 
		pL2 += dim2R; node = brthtab[node];
	}
	dimF = vFact[task].dimF; dim2R = dimF - dimMF; 
	// Define the sectors of the vector in the node task
	ptrO = IdentifyVectorResolution (vFact+node, vidO);
	dimF = vFact[task].dimF; 
	pR2 = ptrO + (vFact[task].dimX-dimF); 
	pR1 = pR2 + dim2R; pR3 = pR1 + dim1;
	// Copy the sectors from the node fthr to the node task
	CopyDoubles (pL1, pR1, dim1 );
	CopyDoubles (pL2, pR2, dim2R);
	CopyDoubles (pL3, pR3, dim3 );
}

/*********************************************************************************/

// This routine computes the vector to be processed, adding the vector vidI computed
// by the children of node task on the vector vidW of node task.
// The vector vidA of node task is used as an auxiliar vector.
void CopyNodeFromChildrenV (ptr_ILU0Factor vFact, int task, int vidI, int vidW, int vidA) {
	// Definition of the global vectors and variables
	int dimM, *chldtab, *nmchtab, *brthtab;
	// Definition of the local vectors and variables
	int node, node1, node2, nmch = 0, even_odd, dimXX, dimM1, dimM2;
	double *pL, *pL1, *pR1, *ptrW, *ptrA;

	// Parameters validation
	if ((vFact == NULL) || (task < 0) || (vidI < 0) || (vidW < 0) || 
			(vidA < 0) || (vFact->mTab == NULL)) {
		printf ("Incorrect parameters in CopyNodeFromChildrenV (%d,%d,%d,%d)\n", 
								task, vidI, vidW, vidA);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	chldtab = vFact[task].mTab[CHLD]; nmchtab = vFact[task].mTab[NMCH]; 
	brthtab = vFact[task].mTab[BRTH]; 
	dimM = vFact[task].dimM; dimXX =  vFact[task].dimX;
	// Create new structures
	if (vFact[task].mVcL == NULL) {
		CreateMatrixDoubles (&(vFact[task].mVcL), SIZE_MAT_VCL, dimXX);
		InitDoubles (*vFact[task].mVcL, SIZE_MAT_VCL*dimXX, -2.0, 0.0);
		vFact[task].dimV = dimXX; vFact[task].tskV = task;
	}
	// Determinate the vectors vidW and vidA of node task 
	ptrW = IdentifyVectorResolution (vFact+task, vidW);
	ptrA = IdentifyVectorResolution (vFact+task, vidA);
	// Acumulate the result matrices of the children
	node = task; node1 = chldtab[node]; node2 = brthtab[node1];
	even_odd = (nmchtab[node] % 2);
	while (node2 != -1) {
		// Define the sizes and the sections of the vectors to be acumulated
		// Definition of the data for the first operand
		if (nmch == 0) {
			dimM1 = vFact[node1].dimM;
			if (chldtab[node1] == -1)
				pR1 = IdentifyVectorResolution (vFact+node1, vidI)+(dimM1-dimM); 
			else
				pR1 = IdentifyVectorResolution (vFact+node1, vidW)+(dimM1-dimM); 
		} else {
			pR1 = (even_odd)? ptrW: ptrA; 
		}
		// Definition of the data for the second operand
		dimM2 = vFact[node2].dimM;
		if (chldtab[node2] == -1)
			pL1 = IdentifyVectorResolution (vFact+node2, vidI)+(dimM2-dimM); 
		else
			pL1 = IdentifyVectorResolution (vFact+node2, vidW)+(dimM2-dimM); 
		// Select where copy the result 
		pL = (even_odd)? ptrA: ptrW; 
		// Acumulate the sectors of the vectors
		CopyDoubles (pR1, pL, dimM); AxpyDoubles (1.0, pL1, pL, dimM); 
		// Prepare the next iteration
		nmch = ((nmch == 0)? 1: nmch) + 1;
		node1 = node2; node2 = brthtab[node1]; even_odd = (1 - even_odd);
	}
}

// This routine computes the vector to be processed, copying the vector vidW
// computed by the father of node task on the vector vidO of node task.
void CopyNodeFromFatherV (ptr_ILU0Factor vFact, int task, int vidW, int vidO) {
	// Definition of the global vectors and variables
	int dimM, *treetab, *chldtab;
	// Definition of the local vectors and variables
	int fthr, dimMF; // , dimXX;
	double *pL1, *pR1;

	// Parameters validation
	if ((vFact == NULL) || (task < 0) || (vidW < 0) || (vidO < 0) || 
			(vFact->mTab == NULL)) {
		printf ("Incorrect parameters in CopyNodeFromFatherV (%d,%d,%d)\n", 
								task, vidW, vidO);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	dimM = vFact[task].dimM; 
	treetab = vFact[task].mTab[TREE]; chldtab = vFact[task].mTab[CHLD];
	// Calculate the properties of the node task
	dimM = vFact[task].dimM; // dimXX = vFact[task].sprM.dim1;
	// Calculate the properties of the father of the node task
	fthr = treetab[task]; dimMF = vFact[fthr].dimM; 
	pL1 = IdentifyVectorResolution (vFact+fthr, vidW); 
	if (chldtab[task] == -1)
		pR1 = IdentifyVectorResolution (vFact+task, vidO) + (dimM-dimMF); 
	else
		pR1 = IdentifyVectorResolution (vFact+task, vidW) + (dimM-dimMF); 
	// Copy the vector from the father of the node task
	CopyDoubles (pL1, pR1, dimMF);
}

/*********************************************************************************/

// This routine computes the vector to be processed, adding the vectors vidI1 computed 
// by the children of the task and the result is included in the vector vidO1 of vFact. 
// In the node task, the vector vidI1 is used as an auxiliar vector.
// It also computes the vector to be processed, adding the vectors vidI2 computed 
// by the children of the task and the result is included in the vector vidO2 of vFact. 
// In the node task, the vector vidA is used as an auxiliar vector.
void BuildCopyNodeFromChildrenV (ptr_ILU0Factor vFact, int task, int vidI1, int vidO1, 
																	int vidI2, int vidO2, int vidA) {
	// Definition of the global vectors and variables
	int dimM, *sizetab, *chldtab, *nmchtab, *brthtab;
	// Definition of the local vectors and variables
	int node, node1, node2, nmch = 0, even_odd;
	int dimXX, dim1, dim2R, dim2L, dim3;
	int dimBF1, dimBF2, dimBM1, dimBM2, dimCM1, dimCM2;
	double *pBL, *pBL1, *pBL2, *pBL3, *pBR1, *pBR2, *pBR3;
	double *pCL, *pCL1, *pCR1, *ptrI1, *ptrO1, *ptrO2, *ptrA;

	// Parameters validation
	if ((vFact == NULL) || (task < 0) || (vidI1 < 0) || (vidO1 < 0) || (vidI2 < 0) || 
			(vidO2 < 0) || (vidA < 0) || (vFact->mTab == NULL)) {
		printf ("Incorrect parameters in BuildCopyNodeFromChildrenV (%d,%d,%d,%d,%d,%d)\n", 
								task, vidI1, vidO1, vidI1, vidO2, vidA);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	dimM = vFact[task].dimM; 
	sizetab = vFact[task].mTab[SIZE]; chldtab = vFact[task].mTab[CHLD];
	nmchtab = vFact[task].mTab[NMCH]; brthtab = vFact[task].mTab[BRTH]; 
	dimXX =  vFact[task].dimX; dim1 = sizetab[task]; dim3 = dimM -dim1;
	// Create new structures
	if ((vFact[task].mVcL == NULL) && (((vidI1 / 10) == 1) || ((vidO1 / 10) == 1) || 
																		 ((vidI2 / 10) == 1) || ((vidO2 / 10) == 1))) {
		CreateMatrixDoubles (&(vFact[task].mVcL), SIZE_MAT_VCL, dimXX);
		InitDoubles (*vFact[task].mVcL, SIZE_MAT_VCL*dimXX, -2.0, 0.0);
		vFact[task].dimV = dimXX; vFact[task].tskV = task;
	}
	if ((vFact[task].mPCG == NULL) && (((vidI1 / 10) == 0) || ((vidO1 / 10) == 0) || 
																		 ((vidI2 / 10) == 0) || ((vidO2 / 10) == 0))) {
		CreateMatrixDoubles (&(vFact[task].mPCG), SIZE_MAT_PCG, dimM);
		InitDoubles (*vFact[task].mPCG, SIZE_MAT_PCG*dimM, -1.0, 0.0);
		vFact[task].dimV = dimXX; vFact[task].tskV = task;
	}
	// Determinate the vectors vidI and vidO of node task 
	ptrI1 = IdentifyVectorResolution (vFact+task, vidI1);
	ptrO1 = IdentifyVectorResolution (vFact+task, vidO1);
	ptrO2 = IdentifyVectorResolution (vFact+task, vidO2);
	ptrA  = IdentifyVectorResolution (vFact+task, vidA );
	// Acumulate the result matrices of the children
	node = task; node1 = chldtab[node]; node2 = brthtab[node1];
	even_odd = (nmchtab[node] % 2);
	while (node2 != -1) {
		// Define the sizes and the sections of the vectors to be acumulated
		// Definition of the data for the first operand
		if (nmch == 0) { // If no accumulation has been made
			dimBM1 = vFact[node1].dimX; dimBF1 = vFact[node1].dimF;
			dim2R = dimBF1 - dimM; 
			pBR2 = IdentifyVectorResolution (vFact+node1, vidI1)+(dimBM1-dimBF1); 
			pBR1 = pBR2 + dim2R; pBR3 = pBR1 + dim1;
		} else { // If some accumulation has been made
			dimBF1 = dimXX; dim2R = dimBF1 - dimM; 
			pBR2 = (even_odd)? ptrO1: ptrI1; 
			pBR1 = pBR2 + dim2R; pBR3 = pBR1 + dim1 ;
		}
		// Definition of the data for the second operand
		dimBM2 = vFact[node2].dimX; dimBF2 = vFact[node2].dimF;
		dim2L = dimBF2 - dimM; 
		pBL2 = IdentifyVectorResolution (vFact+node2, vidI1)+(dimBM2-dimBF2); 
		pBL1 = pBL2 + dim2L; pBL3 = pBL1 + dim1;

		// Select where acumulate the result and its size
		dimXX = dimM+dim2R+dim2L; 
		pBL = (even_odd)? ptrI1: ptrO1; 
		// Acumulate the sectors of the vectors
		CopyDoubles (pBR2, pBL, dim2R); pBL += dim2R;
		CopyDoubles (pBL2, pBL, dim2L); pBL += dim2L;
		CopyDoubles (pBR1, pBL, dim1 ); AxpyDoubles (1.0, pBL1, pBL, dim1); pBL += dim1;
		CopyDoubles (pBR3, pBL, dim3 ); AxpyDoubles (1.0, pBL3, pBL, dim3); pBL += dim3;
#ifdef PRINT_BUILDCOPY_VECTORS
		{
			double *ptr = NULL; int tam = 0;

			sprintf (filename, "BCA_%2.2d_%2.2d", task, node1); 
			ptr = IdentifyVectorResolution (vFact+node1, vidI1)+(dimBM1-dimBF1);
			tam = dimBF1;
			WriteFDoubles (filename, ptr, tam, 40, 30);

			sprintf (filename, "BCA_%2.2d_%2.2d", task, node2); 
			ptr = IdentifyVectorResolution (vFact+node2, vidI1)+(dimBM2-dimBF2);
			tam = dimBF2;
			WriteFDoubles (filename, ptr, tam, 40, 30);

			printf ("(%d) node1 = %d, dimBM1 = %d , dimBF1 = %d , node2 = %d, dimBM2 = %d , dimBF2 = %d\n",
									task, node1, dimBM1 , dimBF1, node2, dimBM2, dimBF2);

			sprintf (filename, "BCA_%2.2d_%2.2d", task, task); 
			ptr = IdentifyVectorResolution (vFact+task , vidO1);
			tam = vFact[task ].dimX;
			WriteFDoubles (filename, ptr, tam, 40, 30);
		}
#endif
		// Define the sizes and the sections of the vectors to be acumulated
		// Definition of the data for the first operand
		if (nmch == 0) {
			dimCM1 = vFact[node1].dimM;
			if (chldtab[node1] == -1)
				pCR1 = IdentifyVectorResolution (vFact+node1, vidI2)+(dimCM1-dimM); 
			else
				pCR1 = IdentifyVectorResolution (vFact+node1, vidO2)+(dimCM1-dimM); 
		} else {
			pCR1 = (even_odd)? ptrO2: ptrA; 
		}
		// Definition of the data for the second operand
		dimCM2 = vFact[node2].dimM;
		if (chldtab[node2] == -1)
			pCL1 = IdentifyVectorResolution (vFact+node2, vidI2)+(dimCM2-dimM); 
		else
			pCL1 = IdentifyVectorResolution (vFact+node2, vidO2)+(dimCM2-dimM); 

		// Select where copy the result 
		pCL = (even_odd)? ptrA: ptrO2; 
		// Acumulate the sectors of the vectors
		CopyDoubles (pCR1, pCL, dimM); AxpyDoubles (1.0, pCL1, pCL, dimM); 
#ifdef PRINT_BUILDCOPY_VECTORS
		{
			double *ptr = NULL; int tam = 0;

			sprintf (filename, "BCB_%2.2d_%2.2d", task, node1); 
			ptr = IdentifyVectorResolution (vFact+node1, vidO2)+(dimCM1-dimM);
			ptr = pCR1;
			tam = dimM;
			WriteFDoubles (filename, ptr, tam, 40, 30);

			sprintf (filename, "BCB_%2.2d_%2.2d", task, node2); 
			ptr = IdentifyVectorResolution (vFact+node2, vidO2)+(dimCM2-dimM);
			ptr = pCL1;
			tam = dimM;
			WriteFDoubles (filename, ptr, tam, 40, 30);

			printf ("(%d) dimCM1 = %d , dimM = %d , dimCM2 = %d , dimM = %d\n",
									task, dimCM1 , dimM, dimCM2, dimM);

			sprintf (filename, "BCB_%2.2d_%2.2d", task, task); 
			ptr = IdentifyVectorResolution (vFact+task , vidO2);
			tam = dimM;
			WriteFDoubles (filename, ptr, tam, 40, 30);
		}
#endif
		// Prepare the next iteration
		nmch = ((nmch == 0)? 1: nmch) + 1;
		node1 = node2; node2 = brthtab[node1]; even_odd = (1 - even_odd);
	}
}

// This routine computes the vector to be processed, copying the vector vidI1
// computed by the father of node task on the vector vidO1 of node task.
// It also computes the vector to be processed, copying the vector vidI2
// computed by the father of node task on the vector vidO2 of node task.
void BuildCopyNodeFromFatherV (ptr_ILU0Factor vFact, int task, int vidI1, int vidO1, 
																int vidI2, int vidO2) {
	// Definition of the global vectors and variables
	int dimM, *treetab, *sizetab, *chldtab;
	int *brthtab = vFact[task].mTab[BRTH];
	// Definition of the local vectors and variables
	int fthr, node, dim1, dim2R, dim3, dimXX, dimF, dimMF;
	double *pL1, *pL2, *pL3, *pR1, *pR2, *pR3, *ptrI, *ptrO;

	// Parameters validation
	if ((vFact == NULL) || (task < 0) || (vidI1 < 0) || (vidO1 < 0) || 
			(vidI2 < 0) || (vidO2 < 0) || (vFact->mTab == NULL)) {
		printf ("Incorrect parameters in BuildCopyNodeFromFatherV (%d,%d,%d,%d,%d)\n", 
								task, vidI1, vidO1, vidI1, vidO2);
		PrintTrace (); exit (-1);
	}
	// Initialization of the vectors
	dimM = vFact[task].dimM; 
	treetab = vFact[task].mTab[TREE]; sizetab = vFact[task].mTab[SIZE];
	chldtab = vFact[task].mTab[CHLD]; brthtab = vFact[task].mTab[BRTH];
	// Calculate the properties of the node
	fthr = treetab[task]; 
	dimMF = vFact[fthr].dimM; dimXX = vFact[fthr].dimX;
	dim1 = sizetab[fthr]; dim3  = vFact[fthr].dimM - dim1; 
	// Define the sectors of the vector in the node fthr
	pL2 = IdentifyVectorResolution (vFact+fthr, vidI1);
	pL3 = pL2 + (dimXX - dim3); pL1 = pL3 - dim1; 
	// Define where the sectors of bad pivots related to node task is
	// and its size
	node = chldtab[fthr]; 
	while (node != task) {
		dimF = vFact[node].dimF; dim2R = dimF - dimMF; 
		pL2 += dim2R; node = brthtab[node];
	}
	dimF = vFact[task].dimF; dim2R = dimF - dimMF; 

	ptrI = IdentifyVectorResolution (vFact+node, vidI1);
	ptrO = IdentifyVectorResolution (vFact+node, vidO1);
	// Define the sectors of the vector in the node task
	dimF = vFact[task].dimF;
	pR2 = ptrO + (vFact[task].dimX-dimF); 
	pR1 = pR2 + dim2R; pR3 = pR1 + dim1;
	// Copy the sectors from the node fthr to the node task
	CopyDoubles (pL1, pR1, dim1 );
	CopyDoubles (pL2, pR2, dim2R);
	CopyDoubles (pL3, pR3, dim3 );

	// Calculate the properties of the father of the node task
	pL1 = IdentifyVectorResolution (vFact+fthr, vidI2); 
	if (chldtab[task] == -1)
		pR1 = IdentifyVectorResolution (vFact+task, vidO2) + (dimM-dimMF); 
	else
		pR1 = IdentifyVectorResolution (vFact+task, vidI2) + (dimM-dimMF); 

	// Copy the vector from the father of the node task
	CopyDoubles (pL1, pR1, dimMF);
#ifdef PRINT_BUILDCOPY_VECTORS
	{
		double *ptr = NULL; int tam = 0;
		char filename[80];

		sprintf (filename, "BCF_%2.2d_%2.2d", fthr, task); 
		ptr = pR1;
		tam = dimMF;
		WriteFDoubles (filename, ptr, tam, 40, 30);

		sprintf (filename, "BCG_%2.2d_%2.2d", fthr, task); 
		ptr = pL1;
		tam = dimMF;
		WriteFDoubles (filename, ptr, tam, 40, 30);

		sprintf (filename, "BCH_%2.2d_%2.2d", fthr, task); 
		if (chldtab[task] == -1)
			ptr = IdentifyVectorResolution (vFact+task, vidO2);
		else
			ptr = IdentifyVectorResolution (vFact+task, vidI2);
		tam = dimM;
		WriteFDoubles (filename, ptr, tam, 40, 30);
	}
#endif
}

/*********************************************************************************/

