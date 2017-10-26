#ifndef ToolsILU0

#define ToolsILU0 1

#include "ILU0Factor.h"

/********************************************************************************************/

extern void FactorILU0SparseMatrix (ptr_ILU0Factor pFact, paramFactor parFac,
                                    int task, int leaf);
/********************************************************************************************/

// This routine applies the preconditioner, going up.
// extern void ILU0ResolutionUp (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf, int task);
extern void ILU0ResolutionUp (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf);
     
// This routine applies the preconditioner, in the top of the tree.
// extern void ILU0ResolutionUpDown (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf, int task);
extern void ILU0ResolutionUpDown (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf);

// This routine applies the preconditioner, going down.
// extern void ILU0ResolutionDown (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf, int task);
extern void ILU0ResolutionDown (ptr_ILU0Factor pFact, int vidI, int vidO, int vidA, int leaf);
/********************************************************************************************/

#endif

