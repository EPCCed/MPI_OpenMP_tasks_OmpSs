#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectorsNew.h"
#include "SparseMatricesNew.h"
#include "SparseHarwellBoeingNew.h"
#include "EliminationTree.h"
#include "TaskQueue.h"
#include "ToolsILU0.h"

#include "Lists.h"
#include "ToolsMPI.h"
#include "ILU0Factor.h"

//#include <blas.h>
//#include <ilupack.h>

//#include <ilupackmacros.h>

/*********************************************************************************/

// This routine computes an elimination tree with nleaves nodes of matrix spr.
// After it creates the vector on which the preconditioner will be,
// filling it with the data related with the tree structure.
// At the end, the filled vector is returned
// The parameter index indicates if 0-indexing or 1-indexing is used.
// The parameter ilpkcomms includes the information related to the communicators 
// on which the computations will be made.
extern ptr_ILU0Factor InitILU0FactorizationMPIOPENMP (SparseMatrix spr, int index, int nleaves,
																								Ilpck_Comm ilpkcomms, char *matrix);

// The routine computes the multilevel ILU factorization of the sparse matrix spr.
// The initial permutation is included in the files whose names appear as parameters.
extern int ILU0FactorizationMPIOPENMP (SparseMatrix spr, int index, int nleaves, paramFactor parFac, 
																						Ilpck_Comm ilpkcomms, ptr_ILU0Factor *ppFact, char *matrix);

