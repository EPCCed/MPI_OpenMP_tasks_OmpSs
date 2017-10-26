#include <stdio.h>
#include <stdlib.h>
#include "InputOutput.h"
#include "SortFunction.h"

/*************************************************************************************/

// Definition of the internal type, useful for this application
typedef struct SortType {
	int *chldtab, *brthtab, *wgthtab, *ownrtab, *hgthtab;
} SortType, *ptr_SortType;

// Creation of the internal type, from the vectors used to define
// the main features of the elimination tree
//   * chldtab[i] includes the first child of the node i
//   * brthtab[i] defines the number of brothers of node i
//   * wgthtab[i] stores the number of nonzeros of the block i
//   * ownrtab[i] indicates the thread on which the node i is mapped
//   * hgthtab[i] informs about the height of the node i
void *CreateSortType (int *chldtab, int *brthtab, int *wgthtab, int *ownrtab, int *hgthtab) {
	ptr_SortType ptr_S;

	// Parameters validation
	if ((chldtab == NULL) || (brthtab == NULL) || (wgthtab == NULL) || 
			(ownrtab == NULL) || (hgthtab == NULL)) {
		printf ("Incorrect parameters in CreateSortType\n"); PrintTrace(); exit (-1);
	}

	// Creation of the local structure
	if ((ptr_S = (ptr_SortType) malloc (sizeof(SortType))) == NULL) {
      printf ("Memory Error in CreateSortType\n"); PrintTrace(); exit (1); 
	}

	// Initialization of the local structure
	ptr_S->chldtab = chldtab; ptr_S->brthtab = brthtab; ptr_S->wgthtab = wgthtab; 
	ptr_S->ownrtab = ownrtab; ptr_S->hgthtab = hgthtab;

	// Return the created structure
	return (void *) ptr_S;
}

// Free the internal dataype included in ptr
void RemoveSortType (void **ptr) {
	ptr_SortType ptr_S;

	if (ptr != NULL) {
		// If a non null pointer has arrived to the routine
		ptr_S = (ptr_SortType) *ptr;
		// If some data can be liberated 
		if (ptr_S != NULL) free (ptr_S); 
		// Initializes correctly the pointer.
		*ptr = NULL;
	}
}

/*************************************************************************************/

// Define if the priority in the node i is greater than in the node j
// Here, it is considered that the priority is defined by levels, 
// having the leaves the greatest priority. The second criterium is 
// the weight of the nodes.
int OrderUpSortType (void *ptr,  int i, int j) {
	ptr_SortType ptrS = (ptr_SortType) ptr;
	int *chldtab, *wgthtab, *hgthtab;
	int res = 0;

	// Parameters validation
	if ((i < 0) || (j < 0) || (ptrS == NULL) || ((chldtab = ptrS->chldtab) == NULL) || 
			((wgthtab = ptrS->wgthtab) == NULL) || ((hgthtab = ptrS->hgthtab) == NULL)) {
		printf ("Incorrect parameters in OrderUpSortType (%d,%d)\n", i, j); PrintTrace(); exit (-1);
	}
	// If both nodes are of the same type
	if (((chldtab[i] == -1) && (chldtab[j] == -1)) ||
			 ((chldtab[i] != -1) && (chldtab[j] != -1))) {
		// To order respect to the height of the node and the weight
		res = ((hgthtab[i] > hgthtab[j]) ||
					 ((hgthtab[i] == hgthtab[j]) && (wgthtab[i] > wgthtab[j])));
	} else 
		// The leaves has greater priority
		res = (chldtab[i] == -1);

	return res;
}

// Define if the priority in the node i is lower than in the node j
// Here, it is considered that the priority is defined by levels, 
// having the leaves the greatest priority. The second criterium is 
// the weight of the nodes.
int OrderDownSortType (void *ptr,  int i, int j) {
	ptr_SortType ptrS = (ptr_SortType) ptr;
	int *chldtab, *wgthtab, *hgthtab;
	int res = 0;

	// Parameters validation
	if ((i < 0) || (j < 0) || (ptrS == NULL) || ((chldtab = ptrS->chldtab) == NULL) || 
			((wgthtab = ptrS->wgthtab) == NULL) || ((hgthtab = ptrS->hgthtab) == NULL)) {
		printf ("Incorrect parameters in OrderDownSortType (%d,%d)\n", i, j); PrintTrace(); exit (-1);
	}
	// If both nodes are of the same type
	if (((chldtab[i] == -1) && (chldtab[j] == -1)) ||
			 ((chldtab[i] != -1) && (chldtab[j] != -1))) {
		// Order respect to the height of the node and the weight
		res = ((hgthtab[i] < hgthtab[j]) ||
					 ((hgthtab[i] == hgthtab[j]) && (wgthtab[i] < wgthtab[j])));
	} else 
		// The leaves has lower priority
		res = (chldtab[j] == -1);

	return res;
}

/*************************************************************************************/

// Seek if the owner of some of the children of the node tsk is tid.
// The function returns true if the node tsk is a leaf, or 
// if the processor tid has solved some children of the node tsk.
int SeekChildren (void *ptr, int tsk, int tid) {
	ptr_SortType ptrS = (ptr_SortType) ptr;
	int *chldtab, *brthtab, *ownrtab;
	int node, res = 1;

	// Parameters validation
	if ((tsk < 0) || (tid < 0) || (ptrS == NULL) || ((chldtab = ptrS->chldtab) == NULL) || 
			((brthtab = ptrS->brthtab) == NULL) || ((ownrtab = ptrS->ownrtab) == NULL)) {
		printf ("Incorrect parameters in SeekChildren (%d,%d)\n", tsk, tid); PrintTrace(); exit (-1);
	}
	node = chldtab[tsk];
	if (node != -1) {
		// If the node has children, the code verifies if the processor tid is 
		// the owner of some child of the node tsk
		while ((node != -1) && (ownrtab[node] != tid)) {
			node = brthtab[node];
		}
		// if the processor tid is not the owner of any children of node tsk,
		// the function results false
		if (node == -1) res = 0;
	}

	return res;
}

// Seek if the owner of the node tsk is tid
int SeekOwner (void *ptr, int tsk, int tid) {
	ptr_SortType ptrS = (ptr_SortType) ptr;
	int *ownrtab;

	// Parameters validation
	if ((tsk < 0) || (tid < 0) || (ptrS == NULL) || ((ownrtab = ptrS->ownrtab) == NULL)) {
		printf ("Incorrect parameters in SeekOwner (%d,%d)\n", tsk, tid); PrintTrace(); exit (-1);
	}
	return (ownrtab[tsk] == tid);
}

/*************************************************************************************/

