#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "reloj.h"
#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectorsNew.h"
#include "SparseMatricesNew.h"
#include "SparseHarwellBoeingNew.h"
#include "EliminationTree.h"
#include "ToolsILU0.h"
#include "ILU0Factor.h"

#include "Lists.h"
#include "ToolsMPI.h"


/*********************************************************************************/

// #define PRINT_TOLERANCE  1
#define ENERGY_NORM  1
#ifdef ENERGY_NORM
#define SIZE_ND_VECT 4
#define EPSILON      1.11E-16
#endif

/*********************************************************************************/

// #define VERBOSE 1

/*********************************************************************************/

// Create the vectors in vFact related to the tasks included in the p_task_queue_goup,
// // where tid is the thread/process related to the queue
extern void InitVectorResolutionMPIOMPSS (ptr_ILU0Factor vFact);

// This routine writes in files the permutation of the leaves.
// The name of the text files are obtained joining fname 
// to the number of the task.
extern void WritePermMPI (ptr_ILU0Factor vFact, char *fname);

// This routine writes in files the vector vid of the leaves.
// The name of the text files are obtained joining fname 
// to the number of the task.
extern void WriteVectorMPIOMPSS (ptr_ILU0Factor vFact, int vid, char *fname);

// This routine reads	in the vector vid of the leaves the files.
// The name of the text files are obtained joining fname 
// to the number of the task.
extern void ReadVectorMPIOMPSS (ptr_ILU0Factor vFact, int vid, char *fname);

// This routine copies the vector vidI to the vector vidO
extern void CopyVectorsMPI (ptr_ILU0Factor vFact, int vidI, int vidO);

extern void InitResolutionMPIOMPSS (ptr_ILU0Factor vFact, double *vec, int vid, int createVect, 
												int nleaves, int adjustVect, Ilpck_Comm ilpkcomms);

void CloseResolutionMPIOMPSS (ptr_ILU0Factor vFact, int nleaves, int vid, double *sol,
                          Ilpck_Comm ilpkcomms);

// Computation of the matrix vector products related to the processor/hebra id,
// on the data included in the vector of preconditioners.
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
extern void ComputeProductMPIOMPSS (ptr_ILU0Factor vFact, int vid1, int vid2);

// #define ERROR_REDUCTION 1
// Computation of the vector operations (optV) related to the processor/hebra id,
// on the data included in the vector of preconditioners.
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
extern double ComputeVectorOperationsMPIOMPSS (ptr_ILU0Factor vFact, int optV, int vid1,
                                    int vid2, double scal);

extern void ReducingOneScalarMPI (double inp, double *out, Ilpck_Comm ilpkcomms);

extern void ReducingTwoScalarsMPI (double inp0, double inp1, double *out0, double *out1,
                              Ilpck_Comm ilpkcomms);

/*********************************************************************************/

// The routine applies the multilevel ILU factorization using the data included in vFact.
// The routine also transforms a consistent to an unconsistent vector.
// Each thread uses a different queue, exploiting the locality.
// The indices vid's specify where the different auxiliar vector is.
extern void ResolutionTransformMPIOMPSS (ptr_ILU0Factor vFact, int nleaves, int vidI, int vidO1, 
																			int vidO2, Ilpck_Comm ilpkcomms, int *val_ordUp, int *tsk_ordUp, int *val_ordDown, int *tsk_ordDown);

/*********************************************************************************/

// The routine computes nItr steps of the PCG on vector rhs, obtaining the vector x,
// by using the preconditioner included in vFact. 
// Initially, the matrix sprR defines the matrix related to the PCG, but it is not used.
// The parameter indexR indicates if 0-indexing or 1-indexing is used.
void ILU0SolverMPIOMPSS (SparseMatrix sprR, int indexR, double *rhs, double *sol, ptr_ILU0Factor vFact, paramFactor parFac, int nleaves, int *itrEnd, double *tolEnd, Ilpck_Comm ilpkcomms);
 
