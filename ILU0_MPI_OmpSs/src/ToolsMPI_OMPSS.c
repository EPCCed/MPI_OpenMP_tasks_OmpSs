#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ScalarVectors.h>
#include <SparseSymmetricNew.h>
#include <reloj.h>
#include <InputOutput.h>
#include "EliminationTree.h"
#include "ToolsMPI.h"
#include "ToolsMPI_OMPSS.h"
#include "Lists.h"


void ILU0LeavesDistribution_OMPSS (SparseMatrix spr, int index, int nleaves, paramFactor parFac,
                                            Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact) {
  // Definition of the global vectors and variables
  int i, j;
  int *v_ord = NULL;
  int remaining_tsks;
  // Definition of the local vectors and variables
  matInts mTab;
  int tsk = -1,  vint[6];
  int my_id, numprocs, root, nleavespr;
  int indexM = 0; 
  MPI_Comm comm;
  double tle1, tle2, tlu1, tlu2, tt1, tt2, tt;
  MPI_Status st;
  int proc, num;
  List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 

  // Initialization of the MPI variables
  srand (0);
  comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
  nleavespr = (nleaves / numprocs);
  if ((nleaves % numprocs) > 0) {
    printf ("Incorrect parameters in ILU0LeavesDistribution (%d,%d)\n", nleaves, numprocs);
		  PrintTrace (); exit (-1);
  }

  // Initialization of auxiliar vectors (mark and ownr are already initialized to 0)
  mTab = vFact->mTab; 
  v_ord = mTab[LEV];

  if (ilpkcomms.color) {
    MPI_Barrier (comm); reloj (&tle1, &tlu1); tt1 = MPI_Wtime ();
    // Initialize the queue of active nodes with the leaves of the tree
    if (my_id == root) { // The root send the data to the other processes
      for (proc = 1; proc < numprocs; proc++) {
				for (j = proc*nleavespr; j < (proc+1)*nleavespr; j++){
					tsk = v_ord[j];
            // Send the data from the original matrix, and remove the arrived packets
            SendLeafFromMatrix (spr, proc, tsk, vFact, &lst, comm);
            CleanList (&lst, TestPacket);
        }
      }
      for (j=0; j<nleavespr; j++) {
        tsk = v_ord[j];
        reloj (&tle1, &tlu1); tt1 = MPI_Wtime ();
        if (vFact[tsk].tLoc == NULL) {
          CreateMatrixDoubles (&vFact[tsk].tLoc, 2, SIZE_TIM_LOC);
          InitDoubles (vFact[tsk].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
        }
        BuildLeafFromMatrixF (spr, indexM, NULL, vFact, tsk);
        reloj (&tle2, &tlu2); tt2 = MPI_Wtime ();
        vFact[tsk].tLoc[0][TLACUM] += 0.0;
        vFact[tsk].tLoc[0][TLBLDP] += (tle2-tle1);
        vFact[tsk].tLoc[1][TLBLDP] += (tlu2-tlu1);
      }
      // Wait until all packets have arrived 
      while (!emptyList (lst)) {
        CleanList (&lst, TestPacket);
      }
    } else {
			for (j = 0; j < nleavespr; j++) {
				tsk = v_ord[my_id*nleavespr + j];
        reloj (&tle1, &tlu1); tt1 = MPI_Wtime ();
        if (vFact[tsk].tLoc == NULL) {
        	CreateMatrixDoubles (&vFact[tsk].tLoc, 2, SIZE_TIM_LOC);
          InitDoubles (vFact[tsk].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
        }
        RecvLeafFromMatrixF (vFact+tsk, root, tsk, comm);
        reloj (&tle2, &tlu2); tt2 = MPI_Wtime ();
        vFact[tsk].tLoc[0][TLACUM] += 0.0;
        vFact[tsk].tLoc[0][TLBLDP] += (tle2-tle1);
        vFact[tsk].tLoc[1][TLBLDP] += (tlu2-tlu1);
      }
    }
  }
}


void ILU0LeavesDistributionPCG_OMPSS (SparseMatrix sprR, int index, int nleaves, paramFactor parFac, Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact) {

  // Definition of the global vectors and variables
  int i,j;
  int *v_ord = NULL;
  int remaining_tsks;
  // Definition of the local vectors and variables
  matInts mTab;
  int tsk = -1, vint[6];
  int my_id, numprocs, root, nleavespr;
  int indexR = 0; 
  MPI_Comm comm;
  double tge1, tge2, tgu1, tgu2, tt1, tt2, tt;
	MPI_Status st;
  int proc, num;
  List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 

  // Initialization of the MPI variables
  srand (0);
  comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
  nleavespr = (nleaves / numprocs);
  if ((nleaves % numprocs) > 0) {
    printf ("Incorrect parameters in ILU0LeavesDistributionPCG (%d,%d)\n", nleaves, numprocs);
    PrintTrace (); exit (-1);
  }
	// Initialization of auxiliar vectors (mark and ownr are already initialized to 0)
  mTab = vFact->mTab; 
  v_ord = mTab[LEV];

  if (ilpkcomms.color) {
    MPI_Barrier (comm); reloj (&tge1, &tgu1); tt1 = MPI_Wtime ();
    // Initialize the queue of active nodes with the leaves of the tree
    if (my_id == root) { // The root send the data to the other processes
      for (proc = 1; proc < numprocs; proc++) {
				for(j = proc*nleavespr; j < (proc+1)*nleavespr; j++){
						tsk = v_ord[j];
						// Send the data from the original matrix, and remove the arrived packets
						SendLeafFromMatrix (sprR, proc, tsk, vFact, &lst, comm);
            CleanList (&lst, TestPacket);
        }
      }
    	for (j = 0; j < nleavespr; j++) {
				tsk = v_ord[j];
        BuildLeafFromMatrixPCG (sprR, indexR, vFact, tsk);
			}
      // Wait until all packets have arrived 
      while (!emptyList (lst)) {
        CleanList (&lst, TestPacket);
      }
    } else {
				for (j = 0; j < nleavespr; j++){
							tsk = v_ord[my_id*nleavespr + j];
							RecvLeafFromMatrixPCG (vFact+tsk, root, tsk, comm);
         }
    }
    MPI_Barrier (comm);
  }
}

// This routine computes the vector operation, defined by optV, which are 
// related to the processor/hebra tid, on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
double ComputeVectorOperationsOMPSS (ptr_ILU0Factor vFact,  
																	int optV, int vid1, int vid2, double scal) {
	// Definition of the local variables
	int tsk, dim, TVect, i, prio;
	double res_gbl = 0.0;
	double *ptr1 = NULL, *ptr2 = NULL;
	double te1, tu1, te2, tu2;
	int dimL =vFact->dimL;
	int *chldtab = vFact->mTab[CHLD];
	int *ownrtab = vFact->mTab[OWNR];

	// Parameters validation
	if ((vFact == NULL) || (optV < 0) || (optV > 5)  
		|| (vid1 < 0) || (vid2 < 0)) {
		printf ("Incorrect parameters in ComputeVectorOperations (%d,%d,%d)\n", 
							optV, vid1, vid2);
		PrintTrace (); exit (-1);
	}
	// Initialize the variable TVect
	switch (optV) {
		case 0: TVect = TLVEC0; break;
		case 1: TVect = TLVEC1; break;
		case 2: TVect = TLVEC2; break;
		case 3: TVect = TLVEC3; break;
		case 4: TVect = TLVEC4; break;
		case 5: TVect = TLVEC5; break;
		default: break;
	}
	  // Create the vectors that will be used to achieve the best order of execution in the leaves

	// Process all the leaves tasks
  for(tsk=0; tsk<dimL; tsk++){
    if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
    #ifdef ORDER
      prio=v_prio[tsk];
    #else
      prio=0;  
    #endif
			#pragma oss task concurrent(res_gbl) priority(prio)  
      {
        double te1, tu1, te2, tu2;
        double res = 0.0;
        reloj (&te1, &tu1);
        // Definition of the vectors required in the computation
        ptr1 = IdentifyVectorResolution (vFact+tsk, vid1);
        ptr2 = IdentifyVectorResolution (vFact+tsk, vid2);
        dim = vFact[tsk].dimM;
			  // Compute the vector operation for the node task        
				switch (optV) {
        case 0: // Copy 
          CopyDoubles  (ptr1, ptr2, dim);      break;
        case 1: // Axpy
          AxpyDoubles (scal, ptr1, ptr2, dim); break;
        case 2: // Dot
          res += DotDoubles (ptr1, ptr2, dim); break;
        case 3: // Xpay
          XpayDoubles (ptr1, scal, ptr2, dim); break;
        case 4: // Vdiv
          VdivDoubles (1.0, ptr1, ptr2, dim);  break;
        case 5: // Vvec
          VvecDoubles (1.0, ptr1, ptr2, 0.0, ptr2, dim); break;
        default: break;
        }
        reloj (&te2, &tu2);
				vFact[tsk].tLoc[0][TVect] += te2-te1; vFact[tsk].tLoc[1][TVect] += tu2-tu1;
				if (optV == 2) {
						#pragma oss critical
          	res_gbl += res;
        }
			}
		}
	}
	#pragma oss taskwait 
	return res_gbl;
}


// This routine computes the products related to the processor/hebra tid,
// on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
void ComputeProductOMPSS (ptr_ILU0Factor vFact,  
											int vid1, int vid2) {
	// Definition of the local variables
	int tsk;
	double *ptr1 = NULL, *ptr2 = NULL;
	double tle1, tlu1, tle2, tlu2; 
	int dimL = vFact->dimL;
  int *chldtab = vFact->mTab[CHLD];
	int *ownrtab = vFact->mTab[OWNR];
  int prio;

	// Parameters validation
	if ((vFact == NULL) || (vid1 < 0) || (vid2 < 0)) {
		printf ("Incorrect parameters in ComputeProduct (%d,%d)\n", vid1, vid2);
		PrintTrace (); exit (-1);
	}
  // Process all the leaves tasks
  for(tsk=0; tsk<dimL; tsk++){
    if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
    #ifdef ORDER
  		prio=v_prio[tsk];
    #else
      prio=0; 
    #endif
    #pragma oss task priority(prio)
      {
        double tle1, tlu1, tle2, tlu2;
        reloj (&tle1, &tlu1);
				// Definition of the vectors required in the computation
				ptr1 = IdentifyVectorResolution (vFact+tsk, vid1);
        ptr2 = IdentifyVectorResolution (vFact+tsk, vid2);
			 	ProdSymSparseMatrixVector3 (vFact[tsk].sprR, vFact[tsk].indR, ptr1, ptr2);
        // Store the computation cost of the routine
        reloj (&tle2, &tlu2);
        vFact[tsk].tLoc[0][TLPROD] += tle2-tle1; vFact[tsk].tLoc[1][TLPROD] += tlu2-tlu1;
      }
    }
  }
  #pragma oss taskwait
}

// This routine copies the local solution on the final vector (sol)
void ILU0CopyToVector_New (ptr_ILU0Factor pFact, int vid, double *sol) {
	// Definition of the local variables
	int *permM = pFact->permM; 
	int dimXX = pFact->dimM, indexM = pFact->indR;
	double *ptr;

	// Parameters validation
	if ((pFact == NULL) || (vid < 0) || (sol == NULL)) {
		printf ("Incorrect parameters in ILU0CopyToVector (%d)\n", vid);
		PrintTrace (); exit (-1);
	}
	// Definition of the vectors required in the operation
	ptr = IdentifyVectorResolution (pFact, vid);
	CopyInvPermuteDoubles (ptr, sol, permM, indexM, dimXX);
}

