#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
//#include <parmetis.h>
#include <ScalarVectors.h>
#include <SparseSymmetricNew.h>
#include <reloj.h>
#include <InputOutput.h>
#include "SortFunction.h"
#include "EliminationTree.h"
#include "ToolsMPI.h"
#include "ToolsMPI_OPENMP.h"
#include "Lists.h"


void ILU0LeavesDistribution_OPENMP (SparseMatrix spr, int index, int nleaves, paramFactor parFac,
                                            Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact) {
  // Definition of the global vectors and variables
  int j;
  int *v_ord = NULL;
  // Definition of the local vectors and variables
  matInts mTab;
  int tsk = -1;
  int my_id, numprocs, root, nleavespr;
  int indexM = 0;
  MPI_Comm comm;
  double tle1, tle2, tlu1, tlu2, tt1, tt2;
  int proc;
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
  v_ord = mTab[LEV];;

  if (ilpkcomms.color) {
    MPI_Barrier (comm); reloj (&tle1, &tlu1); tt1 = MPI_Wtime ();
    // Initialize the queue of active nodes with the leaves of the tree
    if (my_id == root) { // The root send the data to the other processes
      for (proc = 1; proc < numprocs; proc++) {
				for (j = proc*nleavespr; j < (proc+1)*nleavespr; j++){
					tsk = v_ord[j];
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


void ILU0LeavesDistributionPCG_OPENMP (SparseMatrix sprR, int index, int nleaves, paramFactor parFac, Ilpck_Comm ilpkcomms, ptr_ILU0Factor vFact) {

  // Definition of the global vectors and variables
  int j;
  int *v_ord = NULL;
  // Definition of the local vectors and variables
  matInts mTab;
  int tsk = -1;
  int my_id, numprocs, root, nleavespr;
  int indexR = 0; 
  MPI_Comm comm;
  double tge1, tgu1, tt1;
  int proc;
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
double ComputeVectorOperationsOPENMP (ptr_ILU0Factor vFact, int tid, 
																	int optV, int vid1, int vid2, double scal) {
	// Definition of the local variables
	int tsk, TVect;
	double res_gbl = 0.0;
	int dimL =vFact->dimL;
	int *chldtab = vFact->mTab[CHLD];
	int *ownrtab = vFact->mTab[OWNR];

	// Parameters validation
	if ((vFact == NULL) || (optV < 0) || (optV > 5) || 
			(tid < 0) || (vid1 < 0) || (vid2 < 0)) {
		printf ("Incorrect parameters in ComputeVectorOperations (%d,%d,%d,%d)\n", 
							tid, optV, vid1, vid2);
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

	// Process all the leaves tasks
  #pragma omp parallel 
  {
  #pragma omp single
  {
  for(tsk=0; tsk<dimL; tsk++){
      #pragma omp task firstprivate(tsk) shared(res_gbl, vFact)  
      {
    	if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
				int dim;
        double te1, tu1, te2, tu2;
        double res = 0.0;
				double *ptr1 = NULL, *ptr2 = NULL;
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
      	  	res += DotDoubles (ptr1, ptr2, dim); 
        		break;
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
        	  #pragma omp atomic
        	  res_gbl += res;
        }
      }
    }
  }
 }
 } 
	return res_gbl;
}


// This routine computes the products related to the processor/hebra tid,
// on the data included in the vector of preconditioners,
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
void ComputeProductOPENMP (ptr_ILU0Factor vFact, int tid, 
											int vid1, int vid2) {
	// Definition of the local variables
	int tsk;
	int dimL = vFact->dimL;
  int *chldtab = vFact->mTab[CHLD]; 
	int *ownrtab = vFact->mTab[OWNR];

	// Parameters validation
	if ((vFact == NULL) || (tid < 0) || (vid1 < 0) || (vid2 < 0)) {
		printf ("Incorrect parameters in ComputeProduct (%d,%d,%d)\n", tid, vid1, vid2);
		PrintTrace (); exit (-1);
	}
  // Create the vectors that will be used to achieve the best order of execution in the leaves
  // Process all the leaves tasks
  #pragma omp parallel
  {
  #pragma omp single
  {
  for(tsk=0; tsk<dimL; tsk++){
      #pragma omp task firstprivate(tsk) 
      {
    		if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
      		double tle1, tlu1, tle2, tlu2;
					double *ptr1 = NULL, *ptr2 = NULL;
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
	}
	}
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

// The root receives a factor from a slave (src)
// * comm is the communicator in which the messages is received
void RecvFactorFromSlave_New (ptr_ILU0Factor pFact, int src, MPI_Comm comm) {
  int indexF = 0, vint[6];
  PacketNode pcknode;
  MPI_Status st;
  int err;
  // Receive the sizes to the root
	MPI_Recv (vint, 6, MPI_INT, src, Tag_Receive_Dims_Factor_From_Leaf, comm, &st);
  if (vint[1] != -1) { 
    // Creation of a sparse matrix and other structures
    CreateSparseMatrix (&(pFact->sprF), indexF, vint[2], vint[2], vint[3], 0);
    CreateMatrixDoubles (&(pFact->mDia), SIZE_MAT_DIA, vint[4]);
    CreateInts (&(pFact->headL), vint[5]+1);
    CreateDoubles (&(pFact->divs) , vint[5]);
    // Creation of the packet and reception
    MakeFactorStructPacket (*pFact, 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
    err = MPI_Recv (pcknode.ptr, 1, pcknode.pack, src, Tag_Receive_Data_Factor_From_Leaf, comm, &st);
    // Destruction of the packet
		MPI_Type_free (&(pcknode.pack));
    // Adjusting of the structure
		pFact->dimX = pFact->dimF = vint[2]; pFact->nlev = vint[5];
    // TO REPAIR
     pFact->indF = 0; 
  }
}
