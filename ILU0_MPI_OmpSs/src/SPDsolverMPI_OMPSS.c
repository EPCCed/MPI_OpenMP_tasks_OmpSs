#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

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
#include "ToolsMPI_OMPSS.h"
#include "SPDsolverMPI_OMPSS.h"


/*********************************************************************************/

//#define VERBOSE 1

/*********************************************************************************/

// Create the vectors in vFact related to the tasks included in the p_task_queue_goup,
// // where tid is the thread/process related to the queue
void InitVectorResolutionMPIOMPSS (ptr_ILU0Factor vFact) {
	// Definition of the local vectors and variables
  int tsk, dimL = vFact->dimL;
  int *chldtab = vFact->mTab[CHLD];
  int *ownrtab = vFact->mTab[OWNR];

  // Parameters validation
  if (vFact == NULL)  {
	  //printf ("Incorrect parameters in InitVectorResolution (%d)\n", omp_get_thread_num());
	  printf ("Incorrect parameters in InitVectorResolution \n");
    PrintTrace (); exit (-1);
  }
  // Process all the leaves tasks 
  for (tsk=0; tsk<dimL; tsk++){
	  if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
	    #pragma oss task 
      {
	    	// Initialize the scalars
        vFact[tsk].dimV = vFact[tsk].dimM; vFact[tsk].tskV = tsk;
       // Create the PCG vectors to each node task
       CreateMatrixDoubles (&vFact[tsk].mPCG, SIZE_MAT_PCG, vFact[tsk].dimV);
       // Create the local vectors to each node task
       CreateMatrixDoubles (&vFact[tsk].mVcL, SIZE_MAT_VCL, vFact[tsk].dimV);
     }
   }
 }
 #pragma oss taskwait
}


// This routine writes in files the vector vid of the leaves.
// The name of the text files are obtained joining fname 
// to the number of the task.
void WriteVectorMPIOMPSS (ptr_ILU0Factor vFact, int vid, char *fname) {
	// Parameters validation
	if ((vFact == NULL) || (vid < 0) || (fname == NULL)) {
		printf ("Incorrect parameters in WriteVectorMPI (%d)\n", vid);
		PrintTrace (); exit (-1);
	}
		// Definition of the local variables
	int tsk;
	char filename[80];
	double *ptr = NULL;
	int dimL = vFact->dimL;
  int *chldtab = vFact->mTab[CHLD];
  int *ownrtab = vFact->mTab[OWNR];
	for (tsk = 0; tsk < dimL; tsk++) {
		if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
			#pragma oss task
			{
				// Print the corresponding vector
				sprintf (filename, "%s_%2.2d.txt", fname, tsk);
				printf ("Writing %s\n", filename);
				ptr = IdentifyVectorResolution (vFact+tsk, vid);
				WriteFDoubles (filename, ptr, vFact[tsk].dimM, 40, 30);
			}
		}
	}
	#pragma oss taskwait
}

// This routine reads	in the vector vid of the leaves the files.
// The name of the text files are obtained joining fname 
// to the number of the task.
void ReadVectorMPIOMPSS (ptr_ILU0Factor vFact, int vid, char *fname) {
	// Parameters validation
	if ((vFact == NULL) || (vid < 0) || (fname == NULL)) {
		printf ("Incorrect parameters in ReadVectorMPI (%d)\n", vid);
		PrintTrace (); exit (-1);
	}
 		int tsk, dimF;
		char filename[80];
		double *ptr = NULL, *vec = NULL;
		int dimL = vFact->dimL;
  	int *chldtab = vFact->mTab[CHLD];
  	int *ownrtab = vFact->mTab[OWNR];
		for (tsk = 0; tsk < dimL; tsk++) {
			if(chldtab[tsk] == -1 && ownrtab[tsk] >= 0){
				#pragma oss task
				{
					// Print the corresponding vector
					sprintf (filename, "%s_%2.2d.txt", fname, tsk);
					printf ("Reading %s\n", filename);
					ptr = IdentifyVectorResolution (vFact+tsk, vid);
					dimF = ReadDoubles (filename, &vec);
					if (dimF == vFact[tsk].dimM) {
						CopyDoubles (vec, ptr, dimF);
						RemoveDoubles (&vec);
					} else {
						printf ("Incorrect size of file %s in ReadVectorMPI (%d,%d)\n", filename, dimF, vFact[tsk].dimM);
						PrintTrace (); exit (-1);
					}	
				}
			}
		}
		#pragma oss taskwait
}


void InitResolutionMPIOMPSS (ptr_ILU0Factor vFact, double *vec, int vid, int createVect, 
												int nleaves, int adjustVect, Ilpck_Comm ilpkcomms) {
	// Definition of the global vectors and variables
	int i, j;
	int remaining_tsks;
	// Definition of the local vectors and variables
	matInts mTab;
	int tsk = -1,  vint[6];
  int my_id, numprocs, root, nleavespr, salto, stop, start, l, mask, salgo; 
	MPI_Comm comm;
	double tge1, tge2, tgu1, tgu2, tt1, tt2, tt;
//	PacketNode pcknode;
	MPI_Status st;
	int *v_ord = NULL;  
	List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 

	// Initialization of the MPI variables
	srand (0);
	comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
 	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
	nleavespr = (nleaves / numprocs);
	if ((nleaves % numprocs) > 0) {
		printf ("Incorrect parameters in InitResolutionMPI (%d,%d)\n", nleaves, numprocs);
    PrintTrace (); exit (-1);
  }

	// Initialization of auxiliar vectors (mark and ownr are already initialized to 0)
	mTab = vFact->mTab; 
	v_ord = mTab[LEV];

	if (ilpkcomms.color) {
		// Create the local structures for the PCG, if createVect is true
		if (createVect) InitVectorResolutionMPIOMPSS (vFact);
  	// Initialize the queue of active nodes with the leaves of the tree
		if (my_id == root) { // The root send the data to the other processes
			for (i=1; i<numprocs; i++) {
				for(j = i*nleavespr; j < (i+1)*nleavespr; j++){
						tsk = v_ord[j];
						// Send the data from the original matrix, and remove the arrived packets
						SendLeafFromVector3 (vec, i, tsk, vid, vFact+tsk, &lst, comm);
						CleanList (&lst, TestPacket);
					}
			}
			for (j=0; j<nleavespr; j++) {
				tsk = v_ord[j];
				vFact[tsk].indR = 0;
				BuildLeafFromVectorF (vec, vFact, tsk, vid, adjustVect); 
			}
			// Wait until all packets have arrived 
			while (!emptyList (lst)) {
				CleanList (&lst, TestPacket);
			}
		} else {
			for (j=0; j<nleavespr; j++) {
				tsk = v_ord[my_id*nleavespr+j];
				RecvLeafFromVectorF3 (vFact+tsk, root, vid, tsk, adjustVect, comm);
			}
		}
	}
}

void InitStructureMPIOMPSS (ptr_ILU0Factor vFact, int nleaves, Ilpck_Comm ilpkcomms) {
  int my_id, numprocs, nleavespr;
  MPI_Comm comm;
  // Initialization of the MPI variables
  comm = ilpkcomms.comm_prec; 
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
  nleavespr = (nleaves / numprocs);
  if ((nleaves % numprocs) > 0) {
    printf ("Incorrect parameters in InitStructureMPI (%d,%d)\n", nleaves, numprocs);
    PrintTrace (); exit (-1);
  }

  if (ilpkcomms.color && numprocs > 1) {
    int task, node;
    int remaining_tsks;
    int cond, flag, proc, src, dst;
    int *treetab = vFact->mTab[TREE], *hgthtab = vFact->mTab[HGTH];
    int *chldtab = vFact->mTab[CHLD], *brthtab = vFact->mTab[BRTH];
		int *v_ord = vFact->mTab[LEV];
    MPI_Status st;

    remaining_tsks = 2 * nleavespr - 1;
   	task = v_ord[my_id*nleavespr];
		while (remaining_tsks > 1) {
    	task = treetab[task];
      remaining_tsks /= 2;
    }
		proc = my_id; cond = 1; flag = 1;
    while (cond) {
      dst = proc |  flag; src = proc & ~flag;
      node = chldtab[treetab[task]];
      while (node == task) {
 		  	node = brthtab[node];
      }
			MPI_Sendrecv (&vFact[task].dimF, 1, MPI_INT, (proc == dst)?src:dst, 99,
                    &vFact[node].dimF, 1, MPI_INT, (proc == dst)?src:dst, 99, comm, &st);
			task = treetab[task];
      cond = (hgthtab[task] > 1) && ((proc & flag) > 0); flag <<= 1;
    }
  }
}


// This routine copies the vector vid of the leaves to the global vector sol. 
// It also removes the structures created in the initialization of the resolution.
void CloseResolutionMPIOMPSS (ptr_ILU0Factor vFact, int nleaves, int vid, double *sol, 
													Ilpck_Comm ilpkcomms) {
	// Definition of the global vectors and variables
	int i, j, k;
	int *v_ord = NULL;
	int remaining_tsks;
	// Definition of the local vectors and variables
	matInts mTab;
	int tsk = -1, vint[6];
  int my_id, numprocs, root, nleavespr; 
	MPI_Comm comm;
	double tge1, tge2, tgu1, tgu2, tt1, tt2, tt;
	MPI_Status st;
	int proc, num; 

	// Initialization of the MPI variables
	srand (0);
	comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
 	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
	nleavespr = (nleaves / numprocs);
	if ((nleaves % numprocs) > 0) {
		printf ("Incorrect parameters in CloseResolutionMPIOMPSS (%d,%d)\n", nleaves, numprocs);
    PrintTrace (); exit (-1);
  }

	// Initialization of auxiliar vectors (mark and ownr are already initialized to 0)
	mTab = vFact->mTab; 
	v_ord =mTab[LEV];

	if (ilpkcomms.color) {
		if (my_id == root) { // The root send the data to the other processes
			for (proc=1; proc<numprocs; proc++) {
				for (j = proc*nleavespr; j < (proc+1)*nleavespr; j++) {
					tsk = v_ord[j];
					double *ptr = NULL;
					int tam = 0;
					MPI_Status sta; 
					if (vFact[tsk].mVcL == NULL) printf ("(%d) vFact[%d].mVcL == NULL\n", my_id, tsk);
					if (vFact[tsk].mPCG == NULL) printf ("(%d) vFact[%d].mPCG == NULL\n", my_id, tsk);
					ptr = IdentifyVectorResolution (vFact+tsk, vid);
					MPI_Probe (proc, 89, comm, &sta);
					MPI_Get_count (&sta, MPI_DOUBLE, &tam);
					tam = vFact[tsk].dimM;
					MPI_Recv (ptr, tam, MPI_DOUBLE, proc, 89, comm, &sta);
					vFact[tsk].indR = 0; vFact[tsk].indM = 0;
					ILU0CopyToVector_New (vFact+tsk, vid, sol);
				}
			}
			for (j=0; j<nleavespr; j++) {
				tsk = v_ord[j];
				vFact[tsk].indR = 0; vFact[tsk].indM = 0;
				ILU0CopyToVector_New (vFact+tsk, vid, sol);
			}
		} else {
			for (j=0; j<nleavespr; j++) {
				double *ptr = NULL;
				int tam = 0;
				tsk = v_ord[my_id*nleavespr+j];
				ptr = IdentifyVectorResolution (vFact+tsk, vid);
				tam = vFact[tsk].dimM;
				MPI_Send (ptr, tam, MPI_DOUBLE, root, 89, comm);
			}
		}
	}
}


// Computation of the matrix vector products related to the processor/hebra id,
// on the data included in the vector of preconditioners.
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
void ComputeProductMPIOMPSS (ptr_ILU0Factor vFact, int vid1, int vid2) {
	// Parameters validation
	if ((vFact == NULL) || (vid1 < 0) || (vid2 < 0)) {
		printf ("Incorrect parameters in ComputeProductMPIOMPSS (%d,%d)\n", vid1, vid2);
		PrintTrace (); exit (-1);
	}
#ifdef _USE_EXTRAE_
		// Initialization of the event related to matrix vector product
		Extrae_event (EVENTIDENT, EV_PRODUCT);
#endif
		// Compute the product related to tid-th task
		ComputeProductOMPSS (vFact, vid1, vid2);
#ifdef _USE_EXTRAE_
		// Finalization of the event related to matrix vector product
		Extrae_event (EVENTIDENT, EV_NULL);
#endif
}

// #define ERROR_REDUCTION 1
// Computation of the vector operations (optV) related to the processor/hebra id,
// on the data included in the vector of preconditioners.
// The parameters vid1 and vid2 identify the vectors on which the operation will be made.
double ComputeVectorOperationsMPIOMPSS (ptr_ILU0Factor vFact, int optV, int vid1, 
																		int vid2, double scal) {
	double val_gbl = 0.0;
	// Parameters validation
	if ((vFact == NULL) || (optV < 0) || (optV > 5) || 
			(vid1 < 0) || (vid2 < 0)) {
		printf ("Incorrect parameters in ComputeVectorOperationsOMPG (%d,%d,%d)\n", optV, vid1, vid2);
		PrintTrace (); exit (-1);
	}
	{
		double val_loc;
		
#ifdef _USE_EXTRAE_
		// Initialization of the event related to vector operations
		Extrae_event (EVENTIDENT, EV_VECTOPR);
#endif
		// Computation in the corresponding queue
		val_loc = ComputeVectorOperationsOMPSS (vFact, optV, vid1, vid2, scal);
		// Accumulate the local value if a dot operation is being computed
		if (optV == 2) {
				val_gbl += val_loc;
#ifdef ERROR_REDUCTION
			#pragma oss barrier
			{
				if (val_gbl < 0.0) printf ("(%d,%25.15e) --> %25.15e\n", tid, val_loc, val_gbl);
			}
#endif
		}
#ifdef _USE_EXTRAE_
		// Finalization of the event related to vector operations
		Extrae_event (EVENTIDENT, EV_NULL);
#endif
	}
	// Return the value obtained from the vector operation
	return val_gbl;
}

void ReducingOneScalarMPI (double inp, double *out, Ilpck_Comm ilpkcomms) {
MPI_Allreduce (&inp, out, 1, MPI_DOUBLE, MPI_SUM, ilpkcomms.comm_prec);
}

void ReducingTwoScalarsMPI (double inp0, double inp1, double *out0, double *out1, 
															Ilpck_Comm ilpkcomms) {
	double inp[2], out[2];

	inp[0] = inp0 ; inp[1] = inp1;
	MPI_Allreduce (inp, out, 2, MPI_DOUBLE, MPI_SUM, ilpkcomms.comm_prec);
	*out0 = out[0]; *out1 = out[1];
}

/**********************************************************************************/
void SortTsksRes(ptr_ILU0Factor vFact, int tlbld, int tlres, int *val_order, int *tsk_order, Ints_SortPerm_func funcSort){
	int *v_leaf, *tsk_leaf, *v_int, *tsk_int;
	int dimL = vFact->dimL;
	int nleaves=(dimL+1)/2;
	int i;
	double *tBld = NULL, *tRes = NULL;
	int *sorttab = NULL, *hgthtab = NULL, *chldtab = NULL;
	hgthtab = vFact->mTab[HGTH], chldtab = vFact->mTab[CHLD];
	//Write the times of each task int the resolution
	GetLocalTimer (vFact, tlbld, &tBld); GetLocalTimer (vFact, tlres, &tRes);
  AxpyDoubles (1.0, tBld, tRes, vFact->dimL);
  RemoveDoubles (&tBld);
	// Read the vector with the times of each task in the resolutionUp and store in sorttabUp the time fore each task
	if(funcSort == HeapSortPermInts_DecrPos){
		sorttab = vFact->mTab[SORU] ;
	}
	else{
		sorttab = vFact->mTab[SORD] ;
	}
  for (i=0; i<dimL; i++){
		if(chldtab[i] == -1){
	 		sorttab[i] = (-1000000)*tRes[i]; 
		}
		else{
	 		sorttab[i] = 1000000*tRes[i];
		}
	}
  RemoveDoubles (&tRes);
	CreateInts(&v_leaf, dimL); CreateInts(&tsk_leaf, dimL);
	CreateInts(&v_int, dimL); CreateInts(&tsk_int, dimL);
	CopyInts (sorttab, v_leaf, dimL);
  CopyInts (hgthtab, v_int, dimL);
  InitInts(tsk_leaf, dimL, 0, 1);
  InitInts(tsk_int, dimL, 0, 1);
	HeapSortPermInts(v_leaf, tsk_leaf, dimL, HeapSortPermInts_IncrPos);
	HeapSortPermInts(v_int, tsk_int, dimL, funcSort);
	if(funcSort == HeapSortPermInts_DecrPos){
		CopyInts(tsk_leaf, tsk_order, nleaves);
		CopyInts(tsk_int + nleaves, tsk_order + nleaves, nleaves-1);
		CopyInts(v_leaf, val_order, nleaves);
		CopyInts(v_int + nleaves, val_order + nleaves, nleaves-1);
		for(i=0; i<dimL; i++){
			if(chldtab[tsk_order[i]] == -1){
				val_order[i]=val_order[i]*(-1);
			}
		}
	}
	else{
		CopyInts(tsk_int, tsk_order, nleaves-1);
  	CopyInts(tsk_leaf, tsk_order + nleaves-1, nleaves);
		CopyInts(v_int, val_order, nleaves-1);
  	CopyInts(v_leaf, val_order + nleaves-1, nleaves);
		for(i=0; i<dimL; i++){
			if(chldtab[tsk_order[i]] == -1){
				val_order[i]=val_order[i]*(-1);
			}
		}
	}
	RemoveInts(&tsk_leaf); RemoveInts(&tsk_int); RemoveInts(&v_leaf); RemoveInts(&v_int);
}

/**********************************************************************************/
void ILUResUp(ILU0Factor *vFact, int tsk, int tskL, int *chldtab, int **dependencies, int vidI, int vidO2, int vidO1, int vidV, int vidB, int vidS, int vidA, int vidT, int vidX, int *loc_remaining_tsks_godown){
	//Local variables
	double tle1, tlu1, tle2, tlu2;
  int *treetab, *sizetab;
	int chld;
	treetab = vFact->mTab[TREE]; sizetab = vFact->mTab[SIZE];
	// Calculate the properties of the node
	reloj (&tle1, &tlu1);
	chld = chldtab[tsk]; 
	// Obtain the vector on apply the preconditioner
	// Copy the vector to the children
	if (chld != -1) { // The new node is not a leaf
		BuildCopyNodeFromChildrenV (vFact, tsk, vidB, vidV, vidI, vidT, vidX);
	} else if (vidI != vidO2) {
		CopyDoubles  (IdentifyVectorResolution (vFact+tsk, vidI), 
									IdentifyVectorResolution (vFact+tsk, vidO2), 
									sizetab[tsk]);
	}
	reloj (&tle2, &tlu2);
	vFact[tsk].tLoc[0][TLBLDU] += tle2-tle1;
	// Verify if a new task has been obtained
	if (treetab[tsk] != -1) { // the node is not the root
		// Simulation of the solution of the multilevel ILU preconditioner
		reloj (&tle1, &tlu1);
		if (chld == -1) {
			ILU0ResolutionUp (vFact+tsk, vidI, vidB, vidA, (chld == -1));
		} else {
			ILU0ResolutionUp (vFact+tsk, vidV, vidB, vidA, (chld == -1));
		}
		reloj (&tle2, &tlu2);
		vFact[tsk].tLoc[0][TLRESU] += tle2-tle1;
	} else { // The node is the root
		// Simulation of the solution of the multilevel ILU preconditioner
		reloj (&tle1, &tlu1);
		if (chld == -1) {
			ILU0ResolutionUpDown (vFact+tsk, vidI, vidO1, vidA, (chld == -1));
		} else {
			ILU0ResolutionUpDown (vFact+tsk, vidV, vidS, vidA, (chld == -1));
		}
		reloj (&tle2, &tlu2);
		vFact[tsk].tLoc[0][TLRESU] += tle2-tle1;
		*loc_remaining_tsks_godown--;
	}
}

/**********************************************************************************/

#pragma oss task in(chldtab[dependencies[tskL][0]], chldtab[dependencies[tskL][1]]) out(chldtab[tskL]) label(ILURESUp_leaf) priority(prio)
void ILUResUp_leaf(ILU0Factor *vFact, int tsk, int tskL, int *chldtab, int **dependencies, int vidI, int vidO2, int vidO1, int vidV, int vidB, int vidS, int vidA, int vidT, int vidX, int prio, int *loc_remaining_tsks_godown){
  ILUResUp(vFact, tsk, tskL, chldtab, dependencies, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX, loc_remaining_tsks_godown);
}

/**********************************************************************************/

#pragma oss task in(chldtab[dependencies[tskL][0]], chldtab[dependencies[tskL][1]]) out(chldtab[tskL]) label(IluResUp_noleaf) priority(prio)
void ILUResUp_noleaf(ILU0Factor *vFact, int tsk, int tskL, int *chldtab, int **dependencies, int vidI, int vidO2, int vidO1, int vidV, int vidB, int vidS, int vidA, int vidT, int vidX, int prio, int *loc_remaining_tsks_godown){
	ILUResUp(vFact, tsk, tskL, chldtab, dependencies, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX, loc_remaining_tsks_godown);
}

/**********************************************************************************/

void ILUResDown(ILU0Factor *vFact, int tsk, int tskL, int *chldtab, int **dependencies, int father, int vidI, int vidO2, int vidO1, int vidV, int vidB, int vidS, int vidA, int vidT, int vidX){
	//Local variables
	double tle1, tlu1, tle2, tlu2;
	int chld;
	// Obtain the vector on apply the preconditioner
 	// Copy the vector from the father
 	reloj (&tle1, &tlu1);
	BuildCopyNodeFromFatherV (vFact, tsk, vidS, vidB, vidT, vidO2);
	reloj (&tle2, &tlu2);
	vFact[tsk].tLoc[0][TLBLDD] += tle2-tle1;
	// Simulation of the solution of the multilevel ILU preconditioner
	chld = chldtab[tsk];
	reloj (&tle1, &tlu1);
	if (chld == -1) {
		ILU0ResolutionDown (vFact+tsk, vidB, vidO1, vidA, (chld == -1));
	} else {
		ILU0ResolutionDown (vFact+tsk, vidB, vidS, vidA, (chld == -1));
	}
	reloj (&tle2, &tlu2);
	vFact[tsk].tLoc[0][TLRESD] += tle2-tle1;
}

/**********************************************************************************/

#pragma oss task in(chldtab[dependencies[tskL][2]]) out(chldtab[tskL]) label(ILUResDown_leaf) priority(prio)
void ILUResDown_leaf(ILU0Factor *vFact, int tsk, int tskL, int *chldtab, int **dependencies, int father, int vidI, int vidO2, int vidO1, int vidV, int vidB, int vidS, int vidA, int vidT, int vidX, int prio){
	ILUResDown(vFact, tsk, tskL, chldtab, dependencies, father, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX);
}

/**********************************************************************************/

#pragma oss task in(chldtab[dependencies[tskL][2]]) out(chldtab[tskL]) label(ILUResDown_noleaf) priority(prio)
void ILUResDown_noleaf(ILU0Factor *vFact, int tsk, int tskL, int *chldtab, int **dependencies, int father, int vidI, int vidO2, int vidO1, int vidV, int vidB, int vidS, int vidA, int vidT, int vidX, int prio){
	ILUResDown(vFact, tsk, tskL, chldtab, dependencies, father, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX);
}

/**********************************************************************************/

// The routine applies the multilevel ILU factorization using the data included in vFact.
// The routine also transforms a consistent to an unconsistent vector.
// Each thread uses a different queue, exploiting the locality.
// The indices vid's specify where the different auxiliar vector is.
void ResolutionTransformMPIOMPSS (ptr_ILU0Factor vFact, int nleaves, int vidI, int vidO1, int vidO2, 
														 Ilpck_Comm ilpkcomms, int *val_ordUp, int *tsk_ordUp, int *val_ordDown, int *tsk_ordDown) {
	MPI_Comm comm;
	double tge1, tge2, tgu1, tgu2, tt1, tt2; 
	int numprocs, my_id, nleavespr;
	MPI_Status st;
	int j;
	ptr_ILU0Factor vFactL;
	// Initialization of the MPI variables
	srand (0);
	comm = ilpkcomms.comm_prec; 
 	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
	nleavespr = (nleaves / numprocs);
	if ((nleaves % numprocs) > 0) {
		printf ("Incorrect parameters in ResolutionTransformMPIOMPSS (%d,%d)\n", nleaves, numprocs);
    PrintTrace (); exit (-1);
  }
	
	// Parameters validation
	if ((vFact == NULL) || (vFact->mTab == NULL) || 
			(vFact->dimT < 1) || (vidI < 0) || (vidO1 < 0) || (vidO2 < 0)) {
		printf ("Incorrect parameters in ResolutionTransformMPIOMPSS (%d,%d,%d)\n", 
						vidI, vidO1, vidO2); 
		PrintTrace (); exit (-1);
	}
	
	// Initialize the marktab vector
	InitInts (vFact->mTab[MARK], vFact->dimT, 0, 0);
	if (ilpkcomms.color) {
		// Definition of the global vectors and variables
		int dimL, hMPI=0, shift;
		int *treetab, *sizetab, *chldtab, *v_ord;
		int *brthtab, *ownrtab = NULL;
		int *sorttabUp, sorttabDown; 
		int vidV = VECT_VECL, vidB = VECT_BUFL, vidS = VECT_SOLL;
		int vidA = VECT_AUXL, vidT = VECT_TRNF, vidX = VECT_AUXL1;
		// Definition of the local vectors and variables
		int loc_remaining_tsks_goup, loc_remaining_tsks_godown;
		int tsk = -1, chld = -1, i;
		double tle1, tlu1, tle2, tlu2;
	
		int *hgthtab = vFact->mTab[HGTH]; 
		int dst, src, flag = 1, cond, proc, node;
		int OMP_MPI = 1;
		
		int **dependencies, child, n_deps=3, father, task, prio; 

		List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 
		
		loc_remaining_tsks_goup = loc_remaining_tsks_godown = 2 * nleavespr - 1;
		
		// Initialization of the vectors
		dimL = vFact->dimL; treetab = vFact->mTab[TREE]; sizetab = vFact->mTab[SIZE];
		chldtab = vFact->mTab[CHLD]; 
		brthtab = vFact->mTab[BRTH]; ownrtab = vFact->mTab[OWNR];
		v_ord = vFact->mTab[LEV];
		dependencies = (int**) malloc (loc_remaining_tsks_goup * sizeof(int*));
		for (tsk=0; tsk<loc_remaining_tsks_goup; tsk++){
			dependencies[tsk] = (int*) malloc (n_deps * sizeof(int));
			child=chldtab[tsk];
			for(int i=0; i < n_deps-1; i++){
				dependencies[tsk][i]=(child==-1)?tsk:child;
			  if(child!=-1)
				 	child=brthtab[child];
			}
			father = treetab[tsk];
			dependencies[tsk][n_deps-1]=(father==-1)?tsk:father; 
		}
		// To priorize tasks 
		if(tsk_ordUp == NULL){
			CreateInts(&val_ordUp, dimL);
			InitInts(val_ordUp, dimL, 0, 0);
			CreateInts(&tsk_ordUp, dimL);
			InitInts(tsk_ordUp, dimL, 0, 1);
		}
		if(tsk_ordDown == NULL){
			CreateInts(&val_ordDown, dimL);
			InitInts(val_ordDown, dimL, 0, 0);
			CreateInts(&tsk_ordDown, dimL);
			InitInts(tsk_ordDown, dimL, dimL-1, -1);
		}
		shift = v_ord[my_id*nleavespr];
		for (i =0; i<loc_remaining_tsks_goup; i++) { 
			tsk = i + shift;
			if (chldtab[tsk] == -1){
				ILUResUp_leaf(vFact, tsk, i, chldtab, dependencies, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX, prio, &loc_remaining_tsks_godown);
			}
			else{
				ILUResUp_noleaf(vFact, tsk, i, chldtab, dependencies, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX, prio, &loc_remaining_tsks_godown);
			}
		}
		#pragma oss taskwait

		// BEGIN MPI Step
		OMP_MPI = 0; loc_remaining_tsks_goup = hMPI = hgthtab[tsk] - 1;   

		if (hMPI > 0) {
      flag = 1; cond = 1; proc = my_id;
      if ((loc_remaining_tsks_goup > 0) && (OMP_MPI == 0) && (treetab[tsk] != -1)) {
        dst = proc |  flag; src = proc & ~flag;
        cond = (proc & flag) > 0; flag <<= 1;
        if (proc == src) {
          SendVectorResolutionTransform (vFact+tsk, dst, tsk, 1, vidB, vidI, vidT, &lst, comm);
          cond = 0; loc_remaining_tsks_goup = 0; OMP_MPI = 1;
        } else {
         	int node = chldtab[treetab[tsk]];
          while ((node != -1) && (node == tsk)) node = brthtab[node];		// The receiver node is the brother of the sender (tsk)
          vFact[node].dimM = vFact[tsk].dimM-sizetab[tsk]+sizetab[node];
          RecvVectorResolutionTransform (vFact+node, src, node, vidB, vidI, vidT, &lst, comm);
        }
      }
      while (loc_remaining_tsks_goup > 0) {
				tsk = treetab[tsk]; 		// La tarea a procesar es el padre de la que ha enviado la informacion
      	// Calculate the properties of the node
      	reloj (&tle1, &tlu1);
      	chld = chldtab[tsk];  ownrtab[tsk] = 1;
        // Obtain the vector on apply the preconditioner
        // Copy the vector to the children
        if (chld != -1) { // The new node is not a leaf
          BuildCopyNodeFromChildrenV (vFact, tsk, vidB, vidV, vidI, vidT, vidX);
        } else if (vidI != vidO2) {
          CopyDoubles  (IdentifyVectorResolution (vFact+tsk, vidI),
                        IdentifyVectorResolution (vFact+tsk, vidO2),
                        sizetab[tsk]);
        }
        reloj (&tle2, &tlu2);
        vFact[tsk].tLoc[0][TLBLDU] += tle2-tle1;
        // Verify if a new task has been obtained
        if (treetab[tsk] != -1) { // the node is not the root
          // Simulation of the solution of the multilevel ILU preconditioner
          reloj (&tle1, &tlu1);
          if (chld == -1) {
            ILU0ResolutionUp (vFact+tsk, vidI, vidB, vidA, (chld == -1));
          } else {
            ILU0ResolutionUp (vFact+tsk, vidV, vidB, vidA, (chld == -1));
          }
          reloj (&tle2, &tlu2);
          vFact[tsk].tLoc[0][TLRESU] += tle2-tle1;
        } else { // The node is the root
          // Simulation of the solution of the multilevel ILU preconditioner
          reloj (&tle1, &tlu1);
          if (chld == -1) {
            ILU0ResolutionUpDown (vFact+tsk, vidI, vidO1, vidA, (chld == -1));
          } else {
            ILU0ResolutionUpDown (vFact+tsk, vidV, vidS, vidA, (chld == -1));
          }
          reloj (&tle2, &tlu2);
          vFact[tsk].tLoc[0][TLRESU] += tle2-tle1;
          // Adjust the number of processed nodes in the godown queue
          loc_remaining_tsks_godown--;
          node = chldtab[tsk];
					while ( (node != -1) && (ownrtab[node] >= 0) ) node = brthtab[node];
          SendVectorResolutionTransform (vFact+tsk, src, tsk, 0, vidS, vidT, vidT, &lst, comm);
				
          loc_remaining_tsks_godown = 2 * nleavespr - 1;
          cond = (proc & flag) > 0; flag >>= 1;
          dst = proc |  flag; src = proc & ~flag;
        }
        // Adjust the number of processed nodes in the goup queue
		    loc_remaining_tsks_goup--;
				
				if ((loc_remaining_tsks_goup > 0) && (OMP_MPI == 0) && (treetab[tsk] != -1)) {
        	dst = proc |  flag; src = proc & ~flag;
          cond = (proc & flag) > 0; flag <<= 1;
          if (proc == src) {
            SendVectorResolutionTransform (vFact+tsk, dst, tsk, 1, vidB, vidI, vidT, &lst, comm);
          	cond = 0; loc_remaining_tsks_goup = 0; OMP_MPI = 1;
          } else {
           	int node = chldtab[treetab[tsk]];
            while ((node != -1) && (node == tsk)) node = brthtab[node];
            vFact[node].dimM = vFact[tsk].dimM-sizetab[tsk]+sizetab[node];
            RecvVectorResolutionTransform (vFact+node, src, node, vidB, vidI, vidT, &lst, comm);
          }
        }
      }
			loc_remaining_tsks_godown = hMPI;

      if (OMP_MPI) {
        OMP_MPI = 0;
        node = treetab[tsk];
        vFact[node].dimM = vFact[tsk].dimM-sizetab[tsk];
        RecvVectorResolutionTransform (vFact+node, dst, node, vidS, vidT, vidT, &lst, comm);
				
        loc_remaining_tsks_godown = 2 * nleavespr - 1;
        cond = (proc & flag) > 0; flag >>= 1;
        dst = proc |  flag; src = proc & ~flag;
        tsk = treetab[tsk];
      } else
        node = tsk;
			while (hMPI > (hgthtab[tsk])) {
				tsk=chldtab[node];
				while ((tsk != -1) && (ownrtab[tsk] < 0)) tsk = brthtab[tsk]; 
        reloj (&tle1, &tlu1);
        BuildCopyNodeFromFatherV (vFact, tsk, vidS, vidB, vidT, vidO2);
        reloj (&tle2, &tlu2);
        vFact[tsk].tLoc[0][TLBLDD] += tle2-tle1;
        // Simulation of the solution of the multilevel ILU preconditioner
        chld = chldtab[tsk];
        reloj (&tle1, &tlu1);
        if (chld == -1) {
          ILU0ResolutionDown (vFact+tsk, vidB, vidO1, vidA, (chld == -1));
        } else {
          ILU0ResolutionDown (vFact+tsk, vidB, vidS, vidA, (chld == -1));
        }
        reloj (&tle2, &tlu2);
        vFact[tsk].tLoc[0][TLRESD] += tle2-tle1;
        // Adjust the number of processed nodes in the godown queue
       	loc_remaining_tsks_godown--;
        cond = (proc & flag) > 0; flag >>= 1;
        dst = proc |  flag; src = proc & ~flag;
        SendVectorResolutionTransform (vFact+tsk, src, tsk, 0, vidS, vidT, vidT, &lst, comm);
        loc_remaining_tsks_godown = 2 * nleavespr - 1;
				node = tsk; 
      }
    }

		// END MPI Step

		loc_remaining_tsks_godown = 2 * nleavespr - 1;
	
		// Loop if there is some task to be processed
		shift = v_ord[my_id*nleavespr];
		for(i = loc_remaining_tsks_godown-1; i >= 0; i--){
    	tsk = i + shift;
    	father=(treetab[tsk]==-1)?tsk:treetab[tsk];
    	if(father != tsk){  //Process all the task except the root
				if(chldtab[tsk] == -1){
        	prio=0;
        	ILUResDown_leaf(vFact, tsk, i, chldtab, dependencies, father, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX, prio);
      	}
      	else{
        	prio=1;
        	ILUResDown_noleaf(vFact, tsk, i, chldtab, dependencies, father, vidI, vidO2, vidO1, vidV, vidB, vidS, vidA, vidT, vidX, prio);
      	}	
    	}
		}
  	#pragma oss taskwait
  	loc_remaining_tsks_godown=0;
	}
}

/*********************************************************************************/

// The routine computes nItr steps of the PCG on vector rhs, obtaining the vector x,
// by using the preconditioner included in vFact. 
// Initially, the matrix sprR defines the matrix related to the PCG, but it is not used.
// The parameter indexR indicates if 0-indexing or 1-indexing is used.
void ILU0SolverMPIOMPSS (SparseMatrix sprR, int indexR, double *rhs, double *sol, ptr_ILU0Factor vFact, paramFactor parFac, int nleaves, int *itrEnd, double *tolEnd, Ilpck_Comm ilpkcomms) {

	MPI_Comm comm;
  int root, numprocs, my_id;
  MPI_Status st;
  // Initialization of the MPI variables
  srand (0);
  comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
  if ((nleaves % numprocs) > 0) {
 		printf ("Incorrect parameters in ResolutionTransformMPI (%d,%d)\n", nleaves, numprocs);
  	PrintTrace (); exit (-1);
  }
 
	// Definition of the local vectors and variables
	int nItr, itrG;
	double tge1, tge2, tgu1, tgu2, tolG;
#ifdef USE_AFFINITY
 	int *cpuArray = NULL;
#endif

	// Parameters validation
	if ((my_id == 0) && ((sprR.dim1 < 1) || (sprR.vptr == NULL) || (sprR.vpos == NULL) ||
                       (sprR.vval == NULL) || (rhs == NULL) || (sol == NULL) || (vFact == NULL) ||
                       ((indexR & (~1)) != 0))) {
  	printf ("Incorrect parameters in ILU0SolverMPI (%d)\n", indexR);
    PrintTrace (); exit (-1);
  }

	// Fix the number of threads and the number of iterations
	nItr = parFac.maxit;
	// Initialization of the timers
	reloj (&tge1, &tgu1);
	vFact->tGlb[0][TGFPCG] = -tge1; vFact->tGlb[1][TGFPCG] = -tgu1; 
	// Definition of nthreads 
#ifdef USE_AFFINITY
	// Prepare data for mapping
	CreateInts (&cpuArray, nthreads);
	PrepareAffinityDown (cpuArray, nthreads, 8, 4, 2);
#endif

 if (my_id == root) {
    ComputeWeightNodes (sprR, indexR, vFact->perm, vFact->mTab[RANG], vFact->mTab[SIZE],
                          vFact->mTab[WGTH], vFact->dimL);
  }
  MPI_Bcast (vFact->mTab[WGTH], vFact->dimT, MPI_INT, root, comm);

	// Parallel region
	MPI_Barrier(comm); reloj (&tge1, &tgu1);
	{
			// Definition of the local vectors and variables
			int i;
			double betac = 0.0, tolc = 0.0, rhoc = 0.0, alpha = 0.0;
			double betau = 0.0, tolu = 0.0, rhou = 0.0;
#ifdef ENERGY_NORM
			int j = 0, cond = 1;
			double nrm0 = 0.0, sumnu = 0.0, nu = 0.0, *nd = NULL;
			matDoubles data;
#else
			double umbral = 1.0e-6;
#endif

#ifdef _USE_EXTRAE_
			// Initialization of the event related to PCG
			Extrae_event (EVENTIDENT, EV_PCG);
#endif
#ifdef USE_AFFINITY
			// Map the threads if the array cpuArray exists
			if (cpuArray != NULL) MapThreadsOnCores (cpuArray, tid);
#endif
			InitStructureMPIOMPSS (vFact, nleaves, ilpkcomms);
     
      ILU0LeavesDistributionPCG_OMPSS (sprR, indexR, nleaves, parFac, ilpkcomms, vFact);
      
			MPI_Barrier(comm); reloj (&tge1, &tgu1);  

			InitResolutionMPIOMPSS (vFact, rhs, VECT_BU, 1, nleaves, 1, ilpkcomms);			// bu   = rhs
			
			InitResolutionMPIOMPSS (vFact, sol, VECT_XC, 0, nleaves, 0, ilpkcomms);			// xc   = sol
			ComputeProductMPIOMPSS (vFact, VECT_XC, VECT_ZU);														// zu    = Au*xc
#ifdef ENERGY_NORM
			// Create the structures for energy-norm
			CreateMatrixDoubles (&data, SIZE_ND_VECT, nItr+1);
			CreateDoubles (&nd, SIZE_ND_VECT);
			// betau = xc*bu = xc*bu
			betau = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_XC, VECT_BU, 0.0);	  // betau = xc*bu
			// tolu = xc*zu = xc*zu
			tolu = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_XC, VECT_ZU, 0.0);		// tolu = xc*zu
			// tolu = xc*zu = xc*zu
			// [betac,tolc] = [betau,tolu]
			ReducingTwoScalarsMPI (betau, tolu, &betac, &tolc, ilpkcomms);
			if (betac < 0) { printf ("BETAC(norm0) is NULL(%25.10e)\n", betac); }
			if (tolc < 0) { printf ("TOLC(norm0) is NULL(%25.10e)\n", tolc); }
			// nrm0 = 2*betac-tolc
			nrm0 = 2*betac-tolc; sumnu = 0.0; j = 0; cond = 1;
	#ifdef VERBOSE
			printf ("beta0 = %10.5e , tol0 = %10.5e , nrm0 = %10.5e\n", 
											betac, tolc, nrm0);
	#endif
#endif
			ComputeVectorOperationsMPIOMPSS (vFact, 0, VECT_BU, VECT_RU, 0.0);						// ru    = bu
			
			ComputeVectorOperationsMPIOMPSS (vFact, 1, VECT_ZU, VECT_RU, -1.0);					// ru    = ru - zu
			
			ResolutionTransformMPIOMPSS (vFact, nleaves, VECT_RU, VECT_YC, VECT_RC, ilpkcomms, NULL, NULL, NULL, NULL);
			
			ComputeVectorOperationsMPIOMPSS (vFact, 0, VECT_YC, VECT_DC, 0.0);						// dc    = yc
			
			betau = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_RU, VECT_YC, 0.0);		// betau = ru*yc

			tolu  = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_RU, VECT_RC, 0.0);		// tolu  = ru*rc
			
			ReducingTwoScalarsMPI (betau, tolu, &betac, &tolc, ilpkcomms);
			
			if ((betac < 0.0) || (tolc < 0.0)) 
				printf ("Negative Errors after ReducingTwoScalars (%f,%f)\n", betac, tolc);
			tolc  = sqrt (tolc);																										// tolc  = sqrt(tolc)
#ifdef VERBOSE
			printf ("beta0 = %10.5e , tol0 = %10.5e\n", betac, tolc);
#endif

#ifdef ENERGY_NORM
			// Initialize the values of required for the energy norm
			for (i=0; i<(SIZE_ND_VECT-2); i++) data[i][0] = 0.0;
			data[SIZE_ND_VECT-2][0] = betac; data[SIZE_ND_VECT-1][0] = tolc;
#endif
			// Begin the main loop
			i = 0; 
#ifdef ENERGY_NORM
			while ((i < nItr) && (cond)) {
#else
			while ((i < nItr) && (tolc > umbral)) {
#endif
				
				ComputeProductMPIOMPSS (vFact, VECT_DC, VECT_ZU);											    // zu    = Au*dc

				rhou = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_DC, VECT_ZU, 0.0);  // rhou = dc*zu

				ReducingOneScalarMPI (rhou, &rhoc, ilpkcomms);

				if (rhoc < 0) { printf ("RHO is NULL(%d,%25.10e)\n", i, rhoc); }

				rhoc = betac / rhoc;																								  // rhoc = betac / rhoc

#ifdef ENERGY_NORM
				// Adapting the required values for the energy norm
				nu = betac * rhoc; nrm0 += nu; nd[j] = nu; 
				sumnu = AddDoubles (nd, SIZE_ND_VECT);
				if (++j == SIZE_ND_VECT) j = 0; 
#endif

				ComputeVectorOperationsMPIOMPSS (vFact, 1, VECT_DC, VECT_XC, rhoc);			  // xc    = xc + rho * dc
				ComputeVectorOperationsMPIOMPSS (vFact, 1, VECT_ZU, VECT_RU, -rhoc);			  // ru    = ru - rho * zu
				ResolutionTransformMPIOMPSS (vFact, nleaves, VECT_RU, VECT_YC, VECT_RC, ilpkcomms, NULL, NULL, NULL, NULL);
				alpha = betac;																												// alpha = betac

				betau = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_RU, VECT_YC, 0.0); // betac = ru*yc

				tolu  = ComputeVectorOperationsMPIOMPSS (vFact, 2, VECT_RU, VECT_RC, 0.0); // tolc  = ru*rc

				// [betac,tolc] = [betau,tolu]
				ReducingTwoScalarsMPI (betau, tolu, &betac, &tolc, ilpkcomms);

				if ((betac < 0.0) || (tolc < 0.0)) 
					printf ("Negative Errors after ReducingTwoScalars (%f,%f)\n", betac, tolc);

#ifdef ENERGY_NORM
				// Compute the condition related to the energy-norm
				cond = ((i < SIZE_ND_VECT) || (sumnu >= (EPSILON*(nrm0+sumnu))));
#endif
				alpha = betac / alpha; 																								// alpha = betac / alpha

				ComputeVectorOperationsMPIOMPSS (vFact, 3, VECT_YC, VECT_DC, alpha);				// dc    = alpha*dc + yc

				tolc  = sqrt (tolc); 																									// tolc  = sqrt(tolc)
#ifdef PRINT_TOLERANCE
	#ifdef ENERGY_NORM
				printf ("TOLERANCE = %d %20.14e %20.14e %20.14e %20.14e %20.14e\n", i, tolc, nrm0, sumnu, nrm0+sumnu, (EPSILON*(nrm0+sumnu)));
	#else
				printf ("TOLERANCE = %d %20.14e \n", i, tolc);
	#endif
#endif
				i++; 
#ifdef ENERGY_NORM
				// Computation for the energy norm
				data[0][i] = rhoc; data[1][i] = alpha; data[2][i] = betac; data[3][i] = tolc;
#endif
			}
				itrG = i; tolG = tolc;
#ifdef VERBOSE
			printf ("beta1 = %10.5e , tol1 = %10.5e\n", betac, tolc);
			printf ("niters = %d of %d, tol = %10.5e \n", i, nItr, tolc);
#endif
#ifdef ENERGY_NORM
			// Liberation of memory for the energy norm
			RemoveDoubles (&nd); RemoveMatrixDoubles (&data);
#endif
			// Collection of the solution
			CloseResolutionMPIOMPSS (vFact, nleaves, VECT_XC, sol, ilpkcomms);
#ifdef _USE_EXTRAE_
				// Initialization of the event related to PCG
				Extrae_event (EVENTIDENT, EV_NULL);
#endif
		}
	int dim1, nnz1, dim2, nnz2, i;
	GetSizesLeavesPCG (vFact, &dim1, &nnz1);
	
	// Clean the preconditioner structures
	CleanILU0FactorVector (&vFact);
	
	// Free the internal structures
#ifdef USE_AFFINITY
	RemoveInts (&cpuArray);
#endif
	*itrEnd = itrG; *tolEnd = tolG;
	reloj (&tge2, &tgu2);
		vFact->tGlb[0][TGFPCG] += tge2     ; vFact->tGlb[1][TGFPCG] += tgu2     ; 
		vFact->tGlb[0][TGPPCG]  = tge2-tge1; vFact->tGlb[1][TGPPCG]  = tgu2-tgu1; 
}

/*********************************************************************************/

