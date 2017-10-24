//#include "IlupackFactor.h"

void ILU0LeavesDistributionPCG_OPENMP (SparseMatrix sprR, int index, int nleaves, paramFactor parFac, Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact);

void ILU0LeavesDistribution_OPENMP (SparseMatrix spr, int index, int nleaves, paramFactor parFac,
                                            Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact);

// This routine computes the vector operation, defined by optV, which are 
// related to the processor/hebra tid, on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
double ComputeVectorOperationsOPENMP (ptr_ILU0Factor vFact, int tid, 
																	int optV, int vid1, int vid2, double scal);


// This routine computes the products related to the processor/hebra tid,
// on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
void ComputeProductOPENMP (ptr_ILU0Factor vFact, int tid, 
											int vid1, int vid2);

// This routine copies the local solution on the final vector (sol)
void IlupackCopyToVector_New (ptr_ILU0Factor pFact, int vid, double *sol);

void ILU0CopyToVector_New (ptr_ILU0Factor pFact, int vid, double *sol);

// The root receives a factor from a slave (src)
// * comm is the communicator in which the messages is received
void RecvFactorFromSlave_New (ptr_ILU0Factor pFact, int src, MPI_Comm comm);
