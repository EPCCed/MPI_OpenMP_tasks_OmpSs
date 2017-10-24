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
#include "SparseSymmetricNew.h"
#include "SparseHarwellBoeingNew.h"
#include "EliminationTree.h"
#include "TaskQueue.h"
#include "ToolsILU0.h"

#include "Lists.h"
#include "ToolsMPI.h"
#include "ILU0Factor.h"
#include "SPDfactorMPI_OPENMP.h"
#include "SPDsolverMPI_OPENMP.h"

/*********************************************************************************/
int main (int argc, char **argv) {
	int root = 0; 
	double *rhs, *sol;
	SparseMatrix spr;
	paramFactor parFac;
	ptr_ILU0Factor vFact = NULL;

  int my_id, numprocs;

	Ilpck_Comm ilpkcomms;

	int ierr = 0, index = 0, nleaves = 0, itrEnd = 0;
	double tolEnd = 0.0;
	char matrix[100];

  MPI_Init (&argc, &argv);

	// Definition of the variables numprocs and my_id
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs); MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	root = numprocs-1;
	root = 0;
	// Definition of the ILU0 communicators
	CreateILU0Communicator (&ilpkcomms, root, numprocs);

  srand (0);
  // Read the original sparse matrix and verify the parameters
  // Verification the parameters
	if (my_id == root)
  	printf ("Verification of the parameters and Reading of the matrix\n");
  VerificationParameters(argc, argv, &nleaves, &parFac);
	if (my_id == root) {
  	// Read the original sparse matrix
  	CreateSparseMatrixHB2 (argv[NUMBER_PARAMS-1], &spr, 1-index);
  }
	strcpy (matrix, argv[NUMBER_PARAMS-1]);

	if (my_id == root) printf ("Initialization of the PCG vectors\n");
  // Create the vectors
  if (my_id == root) {
  	CreateDoubles (&rhs, spr.dim1); CreateDoubles (&sol, spr.dim1);
    // Initialize the vectors
    InitDoubles (rhs, spr.dim1, 0.0, 0.0); InitDoubles (sol, spr.dim1, 1.0, 0.0);
    ProdSymSparseMatrixVector3 (spr, index, sol, rhs);
  }
  
	// Computation of the preconditioner
	if (my_id == root) printf ("(BEGIN) Computation of the preconditioner\n");
	ierr = ILU0FactorizationMPIOPENMP (spr, index, nleaves, parFac, ilpkcomms, &vFact, matrix);
	MPI_Bcast (&ierr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if ((my_id == root) && (ierr != root)) printf ("ierr = %d\n", ierr);
	if (my_id == root) printf ("(END) Computation of the preconditioner\n");
  
	// Clean the structure of the preconditioner
	CleanILU0FactorVector (&vFact);

	if (ierr != 0) {
    printf ("remaining task is equal to %d\n", ierr);
  } else {

		// Resolution of the linear system
		if(my_id==root) printf ("Computation of the PCG\n");
  	if (my_id == root) InitDoubles (sol, spr.dim1, 0.01, 0.001);
  	srand (my_id);
  	if (my_id == root) InitRandDoubles (sol, spr.dim1, 0.0, 1.0);

		ILU0SolverMPIOPENMP (spr, index, rhs, sol, vFact, parFac, nleaves, &itrEnd, &tolEnd, ilpkcomms);

  // Computation of the error
  	if (my_id == 0) {
  		printf ("(%d,%20.15e) , ", itrEnd, tolEnd);
    	printf ("Error1 = %20.15e , ", spr.dim1 - AddDoubles (sol, spr.dim1) );
    	InitDoubles (rhs, spr.dim1, 1.0, 0.0); AxpyDoubles (-1.0, sol, rhs, spr.dim1);
    	printf ("Error2 = %20.15e , ", sqrt(DotDoubles (rhs, rhs, spr.dim1)));
    	ProdSymSparseMatrixVector3 (spr, index, rhs, sol);
    	printf ("A-Error = %20.15e\n", sqrt(DotDoubles (rhs, sol, spr.dim1)));
  	}

  	{
  	double *tBld = NULL, *tPrc = NULL;
  	GetLocalTimer (vFact, TLBLDP, &tBld); GetLocalTimer (vFact, TLFCTP, &tPrc);
  	AxpyDoubles (1.0, tBld, tPrc, vFact->dimL);
  	if (my_id == 0)
  		MPI_Reduce (MPI_IN_PLACE, tPrc, vFact->dimL, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  	else
  		MPI_Reduce (tPrc, tPrc, vFact->dimL, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  	if (my_id == 0) WriteFDoubles ("TimeFactor.txt", tPrc, vFact->dimL, 40, 30);
  	RemoveDoubles (&tPrc); RemoveDoubles (&tBld);
  	int i;
    matDoubles tGlb, tLoc;
    ptr_ILU0Factor vFactG = NULL;
	
  	if (my_id == 0) {
    	CreateMatrixDoubles (&tGlb, 2, SIZE_TIM_GLB); InitDoubles (tGlb[0], 2*SIZE_TIM_GLB, 0.0, 0.0);
      vFactG = CreateILU0FactorVector (vFact->dimL);
     	for (i=0; i<vFact->dimL; i++) {
     		vFactG[i].tGlb = tGlb;        vFactG[i].mTab = vFact->mTab;
       	vFactG[i].dimT = vFact->dimT; vFactG[i].dimL = vFact->dimL;
       	CreateMatrixDoubles (&(vFactG[i].tLoc), 2, SIZE_TIM_LOC);
       	InitDoubles (vFactG[i].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
     	}
    }
    MPI_Reduce (*(vFact->tGlb), (my_id!=0)?NULL:*tGlb, 2*SIZE_TIM_GLB, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

		CreateMatrixDoubles (&tLoc, 2, SIZE_TIM_LOC); InitDoubles (tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
    for (i=0; i<vFact->dimL; i++) {
    	if (vFact[i].tLoc == NULL)
       	MPI_Reduce (*tLoc, (my_id!=0)?NULL:*(vFactG[i].tLoc), 2*SIZE_TIM_LOC, MPI_DOUBLE, MPI_MAX,
       	             0, MPI_COMM_WORLD);
     	else
     	  MPI_Reduce (*(vFact[i].tLoc), (my_id!=0)?NULL:*(vFactG[i].tLoc), 2*SIZE_TIM_LOC, MPI_DOUBLE,
     	                MPI_MAX, 0, MPI_COMM_WORLD);
    }
    RemoveMatrixDoubles (&tLoc);
 		if (my_id == 0) {
	    PrintTimesILU0Factor (2, vFactG, itrEnd); // WaitChar ();
	 		for (i=0; i<vFact->dimL; i++) {
				vFactG[i].tGlb = NULL; vFactG[i].mTab = NULL;
       	vFactG[i].dimT = vFactG[i].dimL = 0;
     	}
    	RemoveILU0FactorVector (&vFactG);
    	RemoveMatrixDoubles (&tGlb);
 		}
		}
	}
	// Free the the vectors
	if (my_id == root) {
		RemoveDoubles (&sol); RemoveDoubles (&rhs); 
	}

	// Free the preconditioner structures
	RemoveILU0FactorVector (&vFact);

	// Free the sparse matrix 
	if (my_id == root) {
		RemoveSparseMatrix (&spr);
	}

	// Free the communicators
	RemoveILU0Communicator (&ilpkcomms); 

	MPI_Finalize ();
	return 0;
}

