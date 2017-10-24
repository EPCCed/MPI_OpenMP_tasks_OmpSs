#ifndef EliminationTree

#define EliminationTree 1

#include "SparseMatricesNew.h"
#include "ILU0Factor.h"

/*********************************************************************************/

// Compute, recursively, the height of each node of the tree, where
// * chldtab[i] includes the head of the list of the children of node i-th
// * brthtab includes the rest of the lists which are headed in chldtab
// * hgthtab[i] includes the heigth of the node i-th
extern void ComputeHeigthNodes (int *chldtab, int *brthtab, int *hgthtab, int root);

// Compute the sizes which are required to create the structures of the node tsk
// * treetab[i] includes the father of the node i-th 
// * sizetab[i] includes the number of rows in the node i-th 
// * wgthtab[i] includes the number of nonzeros in the node i-th 
// * dimM is the number of rows to be computed
// * nnzM is the number of nonzeros to be computed
// * nlev is the number of nodes to arrive to the root
extern void ComputeSizesNodeEliminationTree (int tsk, int *treetab, int *sizetab, 
																					int *wgthtab, int *dimM, int *nnzM, int *nlev);

// From the treetab, this routine computes other vectors related with the tree
// * chldtab[i] includes the head of the list in which the children of node i-th are
// * brthtab includes the rest of the lists which are headed in chldtab
// * nmchtab[i] includes the size of these lists
// * hgthtab[i] includes the height of each node
// The routine returns the number of nodes of the elimination tree
extern int ComputeEliminationTreeVectors (int *treetab, int *chldtab, int *nmchtab, 
																									int *brthtab, int *hgthtab, int size);

/*********************************************************************************/

// This routine computes the number of nonzeros of each row block of the
// sparse matrix spr. The block structure in defined by the permutation
// perm and the auxiliar vectors rangtab and sizetab, while the result
// is incoporated to the auxiliar vector wgthtab. The size of all the
// auxiliar vectors is dimL.
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void ComputeWeightNodes (SparseMatrix spr, int index, int *perm, int *rangtab, 
																int *sizetab, int *wgthtab, int dimL);

// Scale the differents blocks of the sparse matrix to transform to unconsistent:
// * task define the node in the elimination tree
// * treetab and nmchtab are the required vector related to the elimination tree
// * head is the vector which marks the begin and the end of each level (CSR way)
// * dimN and nlev are the dimensions of the previous vectors
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void AdjustLeafEliminationTree (SparseMatrix spr, int index, int task, int *treetab, 
																int *nmchtab, int dimN, int *head, double *divs, int nlev); 

// Remove the contraint blocks for the leaves which are in the nonminimum block.
// * task define the node in the elimination tree
// * treetab, chldtab and brthtab are the required vector related to the elimination tree
// * head is the vector which marks the begin and the end of each level (CSR way)
// * dimN and nlev are the dimensions of the previous vectors
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void AdjustLeafPCGEliminationTree (SparseMatrix spr, int index, int task, int *treetab,
																int *chldtab, int *brthtab, int dimN, int *head, int nlev);

// From the original data (spr,perm) and the elimination tree information (treetab,rangtab),
// the routine builds the data related to a leaf:
// * sprM is the sparse matrix, and permM is the corresponding permutation.
// * headL is the vector which marks the begin and the end of each level (CSR way)
// The parameters index and indexL indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices spr and sprM.
extern void BuildLeafFromMatrix (SparseMatrix spr, int index, int *perm, int *treetab, 
								int *rangtab, int leaf, SparseMatrix sprM, int indexL, int *permM, int *headL);

// From the matrix spr, the routine computes the leaf related to the node task,
// and the resulted matrix is included in the node task
// The parameter indexM indicates if 0-indexing or 1-indexing is used to extract the factor
// If the parameter tvG is not null, it contains the test vector to be used in 
// the preconditioner computation.
extern void BuildLeafFromMatrixF (SparseMatrix spr, int indexM, double *tvG, 
																	ptr_ILU0Factor vFact, int task);

// From the matrix spr, the routine computes the PCG leaf related to the node task,
// and the resulted matrix is included in the node task.
// The parameter indexR indicates if 0-indexing or 1-indexing is used to extract the factor
extern void BuildLeafFromMatrixPCG (SparseMatrix spr, int indexR, ptr_ILU0Factor vFact, int task);

// This routine computes the matrix to be processed, adding the matrices computed 
// by the children of the task and the result is included in the position task of vFact. 
// The parameter indexM indicates if 0-indexing or 1-indexing is used to accumulate the matrices.
// This routine returns the cost of the addition of the matrices of the children
extern double BuildNodeFromChildren (int indexM, ptr_ILU0Factor vFact, int task);

/*********************************************************************************/

// Scale the differents blocks of the vector to transform to unconsistent:
// * vecL is the vector to be adjusted
// * divs includes the divisors related to each level
// * head is the vector which marks the begin and the end of each level (CSR way)
// * nlev is the dimension of the previous vectors
extern void AdjustVector (double *vecL, double *divs, int *head, int nlev); 

// From the original vector (vec) and the local permutation (permM), 
// the routine builds the local vector (vecL)
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void BuildLeafFromVector (double *vec, int *permM, int index, double *vecL, int dimM);

// From the original vector (vec), the routine extracts the components 
// corresponding to the permutation related to the node task, 
// and it stores them in the vector vid related to the node task.
// The vector will be scaled if the parameter adjustV is nonzero.
extern void BuildLeafFromVectorF (double *vec, ptr_ILU0Factor vFact, int task, 
														int vid, int adjustV);

// This routine computes the vector to be processed, adding the vectors vidI computed 
// by the children of the task and the result is included in the vector vidO of vFact. 
// In the node task, the vector vidI is used as an auxiliar vector.
extern void BuildNodeFromChildrenV (ptr_ILU0Factor vFact, int task, int vidI, int vidO);
/*
// This routine defines the vector on the preconditioner will be applied
// in all the children from the solution (solL) of the task.
// The parameter optP defines the way in which the accumulation is made.
extern void BuildNodeToChildrenV (ptr_ILU0Factor vFact, int task, int optP);
*/
// This routine computes the vector to be processed, copying the vector vidI computed
// by the father of node task on the vector vidO of node task.
extern void BuildNodeFromFatherV (ptr_ILU0Factor vFact, int task, int vidI, int vidO);

/*********************************************************************************/
/*
// From the original vector (vec), 
// the routine computes the local vector where the preconditioner will be applied
extern void CopyLeafFromVectorF (double *vec, ptr_ILU0Factor pFact, int task);
*/
// This routine computes the vector to be processed, adding the vector vidI computed
// by the children of node task on the vector vidW of node task.
// The vector vidA of node task is used as an auxiliar vector.
extern void CopyNodeFromChildrenV (ptr_ILU0Factor vFact, int task, int vidI, int vidW, int vidA);

// This routine computes the vector to be processed, copying the vector vidW
// computed by the father of node task on the vector vidO of node task.
extern void CopyNodeFromFatherV (ptr_ILU0Factor vFact, int task, int vidW, int vidO);

/*********************************************************************************/

// This routine computes the vector to be processed, adding the vectors vidI1 computed 
// by the children of the task and the result is included in the vector vidO1 of vFact. 
// In the node task, the vector vidI1 is used as an auxiliar vector.
// It also computes the vector to be processed, adding the vectors vidI2 computed 
// by the children of the task and the result is included in the vector vidO2 of vFact. 
// In the node task, the vector vidA is used as an auxiliar vector.
extern void BuildCopyNodeFromChildrenV (ptr_ILU0Factor vFact, int task, int vidI1, int vidO1, 
																				int vidI2, int vidO2, int vidA);

// This routine computes the vector to be processed, copying the vector vidI1
// computed by the father of node task on the vector vidO1 of node task.
// It also computes the vector to be processed, copying the vector vidI2
// computed by the father of node task on the vector vidO2 of node task.
extern void BuildCopyNodeFromFatherV (ptr_ILU0Factor vFact, int task, int vidI1, int vidO1, 
																				int vidI2, int vidO2);

#endif

