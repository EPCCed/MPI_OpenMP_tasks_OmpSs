#include "ILU0Factor.h"

void ILU0LeavesDistributionPCG_OMPSS (SparseMatrix sprR, int index, int nleaves, paramFactor parFac, Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact) ;

void ILU0LeavesDistribution_OMPSS (SparseMatrix sprR, int index, int nleaves, paramFactor parFac, Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact) ;

// This routine computes the vector operation, defined by optV, which are 
// related to the processor/hebra tid, on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
double ComputeVectorOperationsOMPSS (ptr_ILU0Factor vFact,  
																	int optV, int vid1, int vid2, double scal);


// This routine computes the products related to the processor/hebra tid,
// on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
void ComputeProductOMPSS (ptr_ILU0Factor vFact,  
											int vid1, int vid2);

// This routine copies the local solution on the final vector (sol)
void ILU0CopyToVector_New (ptr_ILU0Factor pFact, int vid, double *sol);

