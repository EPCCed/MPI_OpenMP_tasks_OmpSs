#ifndef SortFunction

#define SortFunction 1

/*************************************************************************************/

// Creation of the internal type, from the vectors used to define
// the main features of the elimination tree
//   * chldtab[i] includes the first child of the node i
//   * brthtab[i] defines the number of brothers of node i
//   * wgthtab[i] stores the number of nonzeros of the block i
//   * ownrtab[i] indicates the thread on which the node i is mapped
//   * hgthtab[i] informs about the height of the node i
extern void *CreateSortType (int *chldtab, int *brthtab, int *wgthtab, int *ownrtab, int *hgthtab);

// Free the internal dataype included in ptr
extern void RemoveSortType (void **ptr);

/*************************************************************************************/

// Definition of the class of functions SortType_func
typedef int (*SortType_func) (void *, int , int );

// Define if the priority in the node i is greater than in the node j
// Here, it is considered that the priority is defined by levels, 
// having the leaves the greatest priority. The second criterium is 
// the weight of the nodes.
extern int OrderUpSortType (void *ptr,  int i, int j);

// Define if the priority in the node i is lower than in the node j
// Here, it is considered that the priority is defined by levels, 
// having the leaves the greatest priority. The second criterium is 
// the weight of the nodes.
extern int OrderDownSortType (void *ptr,  int i, int j);

/*************************************************************************************/

// Definition of the class of functions Seek_func
typedef int (*Seek_func) (void *, int , int );

// Seek if the owner of some of the children of the node tsk is tid.
// The function returns true if the node tsk is a leaf, or 
// if the processor tid has solved some children of the node tsk.
extern int SeekChildren (void *ptr, int tsk, int tid);

// Seek if the owner of the node tsk is tid
extern int SeekOwner (void *ptr, int tsk, int tid);

/*************************************************************************************/

#endif

