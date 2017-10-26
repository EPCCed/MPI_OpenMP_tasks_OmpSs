#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <blas.h>
#include <ilupack.h>
#include <ilupackmacros.h>
#include <mpi.h>

#include "InputOutput.h"
#include "ScalarVectors.h"
#include "SparseVectorsNew.h"
#include "SparseMatricesNew.h"
#include "SparseHarwellBoeingNew.h"
#include "EliminationTree.h"
#include "TaskQueue.h"
#include "ToolsIlupack.h"

#include "Lists.h"
#include "ToolsMPI.h"
#include "IlupackFactor.h"
#include "SPDfactorMPI.h"

/*********************************************************************************/

/*
#define METIS_CALL 1

// The routine computes the multilevel ILU factorization of the sparse matrix spr.
// The initial permutation is included in the files whose names appear as parameters.
ptr_IlupackFactor IlupackFactorizationMPI (SparseMatrix spr, char *permfile, char *rangfile, 
																						char *treefile, Ilpck_Comm ilpkcomms, int optP) { 
	// Definition of the global vectors and variables
	int i, dim = spr.dim1, dimL; 
	int *perm = NULL, *rangtab = NULL, *treetab = NULL, dimT;
	int *sizetab = NULL, *chldtab = NULL, *nmchtab = NULL, *brthtab = NULL, *wgthtab = NULL;
	int remaining_tsks, *marktab = NULL, *ownrtab = NULL, *hgthtab = NULL;
	IlupackFactor *vFact = NULL;
	// Definition of the local vectors and variables
	matDoubles mTab;
	int task = -1, chld = -1, vint[6];
  int my_id, numprocs, root; 
	MPI_Comm comm;
	double te1, te2, tu1, tu2, tt1, tt2, tt;
	PacketNode pcknode;
	MPI_Status st;

	// Definition of the variables numprocs and my_id
	comm = ilpkcomms.comm_metis; root = ilpkcomms.root_m;
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

#ifdef METIS_CALL
	dimT = 2*numprocs;
	if (my_id == root) CreateInts (&perm, dim); 
//	CreateInts (&rangtab, dimT); CreateInts (&treetab, dimT);
//	InitInts (rangtab, dimT, -1, 0); InitInts (treetab, dimT, -1, 0); 
	CreateMatrixInts (&mTab, SIZE_MAT_TAB, dimT); InitInts (mTab[0], dimT*SIZE_MAT_TAB, 0, 0);
	CreateMatrixDoubles (&tGlb, 2, SIZE_TIM_GLB); InitDoubles (tGlb[0], 2*SIZE_TIM_GLB, 0.0, 0.0);

	rangtab = mTab[RANG]; InitInts (rangtab, dimT, -4, 0);
	treetab = mTab[TREE]; InitInts (treetab, dimT, -6, 0);
	MPI_Barrier (comm); reloj (&te1, &tu1);
	tt = ComputeMetisMPI (spr, root, comm, numprocs, rangtab, treetab, perm);
	reloj (&te2, &tu2);
	printf ("Ordering Time (%2d) = %20.15e = (%20.15e,%20.15e)\n", my_id, tt, te2-te1, tu2-tu1);
	if (my_id == root) {
		printf ("%s --> ", "rangtab"); PrintFInts (rangtab, dimT, 6, 0);
		printf ("%s --> ", "treetab"); PrintFInts (treetab, dimT, 6, 0);
	}
#else	
	// Read perm, rangtab, treetab (as the Metis call)
	if (my_id == root) {
		i    = ReadInts (permfile, &perm); // PrintInts (perm, dimT);
		dimT = ReadInts (rangfile, &rangtab); printf ("%s --> ", rangfile); PrintFInts (rangtab, dimT, 6, 0);
		dimT = ReadInts (treefile, &treetab); printf ("%s --> ", treefile); PrintFInts (treetab, dimT, 6, 0);
	} 
	MPI_Bcast (&dimT, 1, MPI_INT, root, comm);
	if (my_id != root) {
		CreateInts (&rangtab, dimT); CreateInts (&treetab, dimT);
	}
	MPI_Bcast (rangtab, dimT, MPI_INT, root, comm);
	MPI_Bcast (treetab, dimT, MPI_INT, root, comm);
	if (my_id == root) {
		printf ("%s --> ", rangfile); PrintFInts (rangtab, dimT, 6, 0);
		printf ("%s --> ", treefile); PrintFInts (treetab, dimT, 6, 0);
	}
#endif

	if (ilpkcomms.color) {
		comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
  	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

		// Create sizetab, chldtab, nmchtab, brthtab, marktab  vectors
		CreateInts (&sizetab, 8*dimT); InitInts (sizetab, 8*dimT, 0, 0);
		chldtab = sizetab + dimT; nmchtab = chldtab + dimT; brthtab = nmchtab + dimT;
		wgthtab = brthtab + dimT; marktab = wgthtab + dimT; ownrtab = marktab + dimT;
		hgthtab = ownrtab + dimT;
		// Compute sizetab, chldtab, nmchtab, brthtab, hgthtab vectors and the variable dimL (the number of nodes)
		dimL = ComputeEliminationTreeVectors (treetab, chldtab, nmchtab, brthtab, hgthtab, dimT);
		ComputeLengthfromHeader (rangtab, sizetab, dimL); 
		// Compute the maximum number of nonzeros of each node (wgthtab vector)
		if (my_id == root) {
			TransformHeadertoLength (spr.vptr, dim);
			for (i=0; i<dimL; i++) {
				wgthtab[i] = AddPermuteInts (spr.vptr+1, perm+rangtab[i], sizetab[i]);
			}
			TransformLengthtoHeader (spr.vptr, dim);
		}
		MPI_Bcast (wgthtab, dimT, MPI_INT, root, comm);

		srand (my_id);
		MPI_Barrier (comm); reloj (&te1, &tu1); tt1 = MPI_Wtime ();
		// Create and initialize the vector of preconditioners
		vFact = CreateIlupackFactorVector (dimL);
		vFact->perm    = perm   ; vFact->mTab = mTab ;
//		vFact->auxltab = sizetab; vFact->rangtab = rangtab; vFact->treetab = treetab;
		vFact->dimL    = dimL   ; vFact->dimT = dimT ;
		for (i=1; i<dimL; i++) vFact[i] = *vFact;
	
		if (my_id == root) { // The root send the data to the other processes
			// Definition of the local vectors and variables for the root
			int src;
			void *ptr = NULL;
			List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 
			TaskQueue tsk_queue, prc_queue; // Queues required to process the applications of the processors
	
//			printf ("hdrtab   --> "); PrintFSeqInts (0, dimL, 6, 0);
//			printf ("sizetab  --> "); PrintFInts (sizetab, dimL, 6, 0);
	
			// Initialize the queue of active nodes with the leaves of the tree
			CreateTaskQueue (&tsk_queue, dimL, CreateSortType(sizetab,dimT), OrderUpSortType, LookChildren);
			for (i=0; i<dimL; i++)
				if (chldtab[i] == -1) {
					EnQueue(&tsk_queue, i); // PrintQueue (tsk_queue);
				}
			// Initialize the queue of waiting processors 
			CreateTaskQueue (&prc_queue, numprocs, NULL, NULL, NULL);
			// Initialize the number of processes to be processed
			remaining_tsks = dimL;
	
			// Loop if there is some task to be processed
			while (remaining_tsks > 0) {
				// Wait until a message is received from other process
				MPI_Recv (vint, 6, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &st);
				src = st.MPI_SOURCE; task = vint[0];
				switch (st.MPI_TAG) {
					case Tag_Demand_Matrix_From_Root: // if a node is demanded
						// The request is enqueued to be processed after
						EnQueue (&prc_queue, src);
						break;
					case Tag_Receive_Dims_Factor_From_Leaf: // if the process sends a factor
						if (vint[1] != -1) { // If some data are sent, the root have to stored
							// Creation of a sparse matrix and reception of the data on it
							CreateSparseMatrix (&(vFact[task].sprF), vint[2], vint[2], vint[3], 0);
							CreateDoubles (&(vFact[task].diag), 3*vint[4]+vint[5]); 
							vFact[task].diag1 = vFact[task].diag  + vint[4];
							vFact[task].diag2 = vFact[task].diag1 + vint[4];
							vFact[task].divs  = vFact[task].diag2 + vint[4];
//						MakeFactorStructPacket2 (vFact[task], 0, vint[2], vint[3], vint[5], &pcknode);
							MakeFactorStructPacket (vFact[task], 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
							MPI_Recv (pcknode.ptr, 1, pcknode.pack, src, Tag_Receive_Data_Factor_From_Leaf, 
													comm, &st);
							MPI_Type_free (&(pcknode.pack));
						}
						// Verify if the father of the task is prepared to be processed
						IlupackEnQueueUp (vFact+task, &tsk_queue, task);
						// Adjust the number of processed nodes
         		remaining_tsks--;
						break;
					default:
						break;
				}
				// If some processes require data and some data exist
				if (!(EmptyQueue(tsk_queue)) && !(EmptyQueue(prc_queue))) {
					// Visit all the requests, assigning a node tsk or reinserting the request
					// The idea is maintain the order of the requests
					for (i=0; i<NumNodesQueue(prc_queue); i++) {
						// Get a processor and a node task available for this processor
						src = DeQueue (&prc_queue, src); task = DeQueue (&tsk_queue, src);
						if (task == -1) // Enqueue again if no task is founded
							EnQueue (&prc_queue, src);
						else { // A node task for src is available
							// The task id is sent to src
	 						ownrtab[task] = src; chld = chldtab[task];
							MPI_Send (&task, 1, MPI_INT, src, Tag_Send_Task_To_Leaf, comm);
							if (chld == -1) { // The new node is a leaf
								// Send the data from the original matrix, and remove the arrived packets
								SendLeafFromMatrix2 (spr, src, task, vFact, &lst, comm);
								CleanList (&lst, TestPacket);
							} else { // The new node is not a leaf
								// Build the node from the factors computed for the children of node task
								BuildNodeFromChildren (vFact, task, optP);
								vint[0] = vFact[task].dimM; vint[1] = vFact[task].sprL.dim1;
								vint[2] = vFact[task].sprL.vptr[vint[1]];
								vint[3] = vFact[task].dimM-sizetab[task];
								vint[4] = vFact[task].nlev;
								// Send the features of the matrix and send the data
								MPI_Send (vint, 5, MPI_INT, src, Tag_Send_Dims_Matrix_To_Leaf, comm);
//							MakeFactorStructPacket2 (vFact[task], 1, vint[1], vint[2], vint[4], &pcknode);
								MakeFactorStructPacket (vFact[task], 1, vint[1], vint[2], vint[3], vint[4], &pcknode);
								MPI_Send (pcknode.ptr, 1, pcknode.pack, src, Tag_Send_Data_Matrix_To_Leaf, comm);
								MPI_Type_free (&(pcknode.pack));
							}
						}
					}
				}
			}

			// Wait until all packets have arrived 
			while (!emptyList (lst)) {
				CleanList (&lst, TestPacket);
			}
			// Send a message to conclude the computation
			for (i=0; i<numprocs; i++)
				if (i != root) { // The root sends the messages 
					if (!LookTaskInQueue(&prc_queue, i))
						// If the processor is still computing, the root have to wait until this computation ends
						MPI_Recv (vint, 6, MPI_INT, i, Tag_Demand_Matrix_From_Root, comm, &st);
					// Inform that the computation is finalized
					MPI_Send (&task, 1, MPI_INT, i, Tag_End_Distribution_To_Leaf, comm);
				}
			// Remove all the nodes of the queue
			for (i=0; i<NumNodesQueue(prc_queue); i++)
				src = DeQueue (&prc_queue, src);
			// Free the internal structures
			ptr = RemoveTaskQueue (&prc_queue); RemoveSortType (&ptr);
			ptr = RemoveTaskQueue (&tsk_queue); RemoveSortType (&ptr); 
	
		}	else {
			// Demand a task to the root
			MPI_Send (&task, 1, MPI_INT, root, Tag_Demand_Matrix_From_Root, comm);
			// Receive the answer from the root
			MPI_Recv (&task, 1, MPI_INT, root, MPI_ANY_TAG, comm, &st);
			// If the root send a node task, it has to be processed
			while (st.MPI_TAG == Tag_Send_Task_To_Leaf) {
				chld = chldtab[task];
				// Receive the matrix from the root
				if (chld == -1) { // By blocks, if the node task is a leaf
					RecvLeafFromMatrixF (vFact+task, root, task, comm);
				} else { // From a chunck, if the node task is not a leaf
					// First, receive the sizes useful to create the local estructures
					MPI_Recv (vint, 6, MPI_INT, root, Tag_Send_Dims_Matrix_To_Leaf, comm, &st);
					// Update the scalars, vectors and structures, which are required to the task
					// Scalars
					vFact[task].dimM = vint[0]; vFact[task].nlev = vint[4];
					// Sparse matrix
					CreateSparseMatrix (&(vFact[task].sprL), vint[1], vint[1], vint[2], 0);
					// Real vectors
					CreateDoubles (&(vFact[task].diag), 3*vint[3]+vint[4]); 
					vFact[task].diag1 = vFact[task].diag  + vint[3];
					vFact[task].diag2 = vFact[task].diag1 + vint[3];
					vFact[task].divs  = vFact[task].diag2 + vint[3];
					// Integer vectors
					CreateInts (&(vFact[task].permL), vint[4]+1); vFact[task].headL = vFact[task].permL;
					// Creation, reception, and freeing of the sparse matrix
//				MakeFactorStructPacket2 (vFact[task], 1, vint[1], vint[2], vint[4], &pcknode);
					MakeFactorStructPacket (vFact[task], 1, vint[1], vint[2], vint[3], vint[4], &pcknode);
					MPI_Recv (pcknode.ptr, 1, pcknode.pack, root, Tag_Send_Data_Matrix_To_Leaf, comm, &st);
					MPI_Type_free (&(pcknode.pack));
				}
	
				// Simulation of the computation of the multilevel ILU preconditioner
				FactorizeSparseMatrix (vFact+task, task); // WaitDelay (wgthtab[task]);
				// Adjust of the scalar values in the structure
				vFact[task].dimX = vFact[task].sprL.dim1; vFact[task].dimF = vFact[task].sprF.dim1;
	
				// Define the sizes to be sent to the root
				vint[0] = task; vint[1] = -1; vint[2] = -1; vint[3] = -1; vint[4] = -1; vint[5] = -1;
				if (vFact[task].sprF.dim1 > 0) { // If a factor has been computed
					vint[1] = vFact[task].dimM;
					vint[2] = vFact[task].sprF.dim1;
					vint[3] = vFact[task].sprF.vptr[vint[2]];
					vint[4] = vFact[task].dimM-sizetab[task];
					vint[5] = vFact[task].nlev;
				} 
				// Send the sizes to the root
				MPI_Send (vint, 6, MPI_INT, root, Tag_Receive_Dims_Factor_From_Leaf, comm);
				if (vint[1] != -1) { // Send the factor to the root, if it exists
//				MakeFactorStructPacket2 (vFact[task], 0, vint[2], vint[3], vint[5], &pcknode);
					MakeFactorStructPacket (vFact[task], 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
					MPI_Send (pcknode.ptr, 1, pcknode.pack, root, Tag_Receive_Data_Factor_From_Leaf, comm);
					MPI_Type_free (&(pcknode.pack));
				}
				// Demand a task to the root
				MPI_Send (&task, 1, MPI_INT, root, Tag_Demand_Matrix_From_Root, comm);
				// Receive the answer from the root
				MPI_Recv (&task, 1, MPI_INT, root, MPI_ANY_TAG, comm, &st);
			}
		}
		// The vector ownrtab is broacasted to all the processors
		MPI_Bcast (ownrtab, dimT, MPI_INT, root, comm);
		reloj (&te2, &tu2); tt2 = MPI_Wtime ();
//		if (my_id == root) 
			printf ("Preconditioner Time (%2.2d) = %20.15e = (%20.15e,%20.15e)\n", my_id, tt2-tt1, te2-te1, tu2-tu1);
		
		// Print ownrtab
		if (my_id == root) {
			printf ("hdrtab   --> "); PrintFSeqInts (0, dimL, 6, 0);
			printf ("ownrtab  --> "); PrintFInts (ownrtab, dimL, 6, 0);
			printf ("perm     --> "); PrintFInts (perm   , dimL, 6, 0);
		}

	}
	
	// Return data
	return vFact;
}
*/
/*********************************************************************************/
/*
// Initialization of the vector of preconditioners,
// where the root is the processor where the data is
void InitTaskQueuesResolutionMPI (ptr_IlupackFactor vFact, int root, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dimL = vFact->dimL, dimT = vFact->dimT;
	int *sizetab = vFact->auxltab;
	int *chldtab = sizetab + dimT, *nmchtab = chldtab + dimT, *brthtab = nmchtab + dimT; 
	int *wgthtab = brthtab + dimT, *marktab = wgthtab + dimT, *ownrtab = marktab + dimT;
	// Definition of the global vectors and variables
	int i, *ntsks = NULL;
	TaskQueue *vtsk_queue_goup, *vtsk_queue_godown;
 	int my_id;

	// Definition of the variable my_id
 	MPI_Comm_rank(comm, &my_id);

	if (my_id == root) {
		//Create the number of tasks of the thread
		CreateInts (&ntsks, 1); InitInts (ntsks, 1, 0, 0);
		// Initialize the number of leaves of the tree
		vFact->ntsks = ntsks;
		for (i=0; i<dimL; i++) {
//			vFact[i].ntsks  = ntsks;
			if (chldtab[i] == -1) (*ntsks)++;
		}
	} else {
		// Create the local queue of active nodes 
		vtsk_queue_goup   = (ptr_TaskQueue) malloc (sizeof(TaskQueue));
		vtsk_queue_godown = (ptr_TaskQueue) malloc (sizeof(TaskQueue));
		CreateTaskQueue (vtsk_queue_goup  , dimL, CreateSortType(sizetab,dimT), OrderUpSortType, NULL);
		CreateTaskQueue (vtsk_queue_godown, dimL, CreateSortType(sizetab,dimT), OrderDownSortType, NULL);

		//Create the number of tasks of the thread
		CreateInts (&ntsks, 1); InitInts (ntsks, 1, 0, 0);
		// Initialize the queue of active nodes with the leaves of the tree, 
		// and the number of task related to the processor
		vFact->goup = vtsk_queue_goup; vFact->godown = vtsk_queue_godown; vFact->ntsks = ntsks;
		for (i=0; i<dimL; i++) {
//		vFact[i].goup = vtsk_queue_goup; vFact[i].godown = vtsk_queue_godown; vFact[i].ntsks = ntsks;
			if (ownrtab[i] == my_id) {
				if (chldtab[i] == -1) EnQueue(vtsk_queue_goup, i); 
		  	(*ntsks)++;
			}
		}
		// Initialize the reset value for each queue
		SetResetQueue(vtsk_queue_goup);
	} 
} 

// Initialization of the vector of preconditioners from the original vector (vec),
// where the root is the processor where the data is
void InitResolutionMPI (ptr_IlupackFactor vFact, double *vec, int root, MPI_Comm comm) {
	int my_id, task, remaining_tsks_goup, dim = 0;
	double *vecL = NULL, *bufL = NULL, *solL = NULL, *auxL = NULL;
		
	// Definition of the variable my_id
  MPI_Comm_rank(comm, &my_id);

	// If the structures are not initialized
	if (vFact->ntsks == NULL) {
 		InitTaskQueuesResolutionMPI (vFact, root, comm);
	}
	// Initialize the number of tasks to be processed
	remaining_tsks_goup = vFact->ntsks[0];

	if (my_id == root) { // The root sends the vectors to the other processors
		int src;
		MPI_Status st;
		List lst = {NULL, NULL}; // List of nonfinalized MPI_Isend 
		// Loop if there is some task to be processed
		while (remaining_tsks_goup > 0) {
			// Wait until a message is received from other process
			MPI_Recv (&task, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &st);
			src = st.MPI_SOURCE; 
			switch (st.MPI_TAG) {
				case Tag_Demand_Vector_From_Root: // If a vector is required
					// Send the vector, and adjust the number of tasks
					SendLeafFromVector2 (vec, src, task, vFact, &lst, comm);
					CleanList (&lst, TestPacket);
       		remaining_tsks_goup--;
					break;
				default:
					break;
			}
		}
		// Wait until all packets have arrived 
		while (!emptyList (lst)) {
			CleanList (&lst, TestPacket);
		}
	} else {
		TaskQueue *vtsk_queue_goup = vFact->goup; // Queue where the local leaves are

		InitVectorResolution (vFact, vtsk_queue_goup, my_id);

		// Initialize the queue to its reset_value
		ResetQueue (vtsk_queue_goup, -1);
		// Process all the task included in the queue
		task = DeQueue (vtsk_queue_goup, my_id);
		while (task != -1) {
			// Demand and receive the vector related to the node task
			MPI_Send (&task, 1, MPI_INT, root, Tag_Demand_Vector_From_Root, comm);
			RecvLeafFromVectorF (vFact+task, root, task, comm);
			task = DeQueue (vtsk_queue_goup, my_id);
		}
	}
}

// Copy to the result vector (sol) from the vector of preconditioners,
// where the root is the processor where the data is
void CloseResolutionMPI (ptr_IlupackFactor vFact, double *sol, int root, MPI_Comm comm) {
 	int my_id, task, remaining_tsks_godown;
	double *vecL = NULL, *bufL = NULL, *solL = NULL, *auxL = NULL;
		
	// Definition of the variable my_id
  MPI_Comm_rank(comm, &my_id);

	// If the structures are not initialized
	if (vFact->ntsks == NULL) {
 		InitTaskQueuesResolutionMPI (vFact, root, comm);
	}
	// Initialize the number of tasks to be processed
	remaining_tsks_godown = vFact->ntsks[0];

	if (my_id == root) { // The root receives the vectors from the other processors
 		int src;
		MPI_Status st;
		// Loop if there is some task to be processed
		while (remaining_tsks_godown > 0) {
			// Wait until a message is received from other process
			MPI_Recv (&task, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &st);
			src = st.MPI_SOURCE; 
			switch (st.MPI_TAG) {
				case Tag_Send_Task_To_Root: // If a vector has to be received
					MPI_Recv (vFact[task].solL, vFact[task].dimM, MPI_DOUBLE, src, Tag_Send_Solution_To_Root, 
											comm, &st);
					// Copy the received vector to the result vector (sol)
					IlupackCopyToVector (vFact+task, NULL, sol);
       		remaining_tsks_godown--;
					break;
				default:
					break;
			}
		}
	} else {
		TaskQueue *vtsk_queue_goup = vFact->goup; // Queue where the local leaves are
	
		// Initialize the queue to its reset_value
		ResetQueue (vtsk_queue_goup, -1);
		// Process all the task included in the queue
		task = DeQueue (vtsk_queue_goup, my_id);
		while (task != -1) {
			// Send the number and the computed vector related to the node task
			// If the task if a leaf (only 1 node in the tree) the solution is sent to the root
			MPI_Send (&task, 1, MPI_INT, root, Tag_Send_Task_To_Root, comm);
			MPI_Send (vFact[task].solL, vFact[task].dimM, MPI_DOUBLE, root, Tag_Send_Solution_To_Root, 
									comm);
			task = DeQueue (vtsk_queue_goup, my_id);
		}

		CloseVectorResolution (vFact, vtsk_queue_goup, my_id);
	}
}

// Simulate the product on the data included in the vector of preconditioners.
void SimulateProductMPI (ptr_IlupackFactor vFact, int root, MPI_Comm comm) {
	int my_id;
		
	// Definition of the variable my_id
 	MPI_Comm_rank(comm, &my_id);
	// Only useful if the processor is not the root
//	if (my_id != root) {
		SimulateProduct (vFact, vFact->goup, my_id);
//	}
}

// Simulate the vector operations (optV) related to the processor/hebra id,
// on the data included in the vector of preconditioners.
void SimulateVectorOperationsMPI (ptr_IlupackFactor vFact, int root, MPI_Comm comm, int optV, double *scal) {
	int my_id;
	double scal_loc = 0.0;
		
	// Definition of the variable my_id
 	MPI_Comm_rank(comm, &my_id);
	// Only useful if the processor is not the root
//	if (my_id != root) {
		scal_loc = SimulateVectorOperations (vFact, vFact->goup, my_id, optV);
//	}
	// Compute the global reduction if it is necessary
	if (optV == 3) {
//		if (my_id != root) printf ("scal_loc[%d] = %f\n", my_id, scal_loc);
		MPI_Allreduce (&scal_loc, scal, 1, MPI_DOUBLE, MPI_SUM, comm);
//		if (my_id == root) printf ("scal_loc[%d] = %f\n", my_id, *scal);
	}
}

// The routine applies the multilevel ILU factorization using the data included in vFact.
// The parameter optP defines the way in which the accumulation is made.
void IlupackResolutionMPI (ptr_IlupackFactor vFact, int root, MPI_Comm comm, int optP) { 
	// Definition of the global vectors and variables
	int dimL = vFact->dimL, dimT = vFact->dimT;
	int *perm = vFact->perm, *treetab = vFact->treetab; 
	int *sizetab = vFact->auxltab;
	int *chldtab = sizetab + dimT, *nmchtab = chldtab + dimT, *brthtab = nmchtab + dimT; 
	int *wgthtab = brthtab + dimT, *marktab = wgthtab + dimT, *ownrtab = marktab + dimT;
	// Definition of the local vectors and variables
	int i, remaining_tsks = 0, remaining_tsks_goup, remaining_tsks_godown;
  int my_id, task = -1, chld = -1, flag, src, dim, vint[6];
	double *vecL = NULL, *bufL = NULL, *solL = NULL, *auxL = NULL;
	TaskQueue *vtsk_queue_goup, *vtsk_queue_godown;
	ptr_SimpleNode smpnode;
	List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 
	MPI_Status st;

	// Definition of the variable my_id
  MPI_Comm_rank(comm, &my_id);

	// Only useful if the processor is not the root
//	if (my_id != root) {
		// Initialize the marktab vector
		InitInts (marktab, dimT, 0, 0);
		// Use the inicialize queue
		vtsk_queue_goup   = vFact->goup; vtsk_queue_godown = vFact->godown;
		// Use the initial value for each queue
		ResetQueue(vtsk_queue_goup, -1);
		// Initialize the number of processes to be processed
		remaining_tsks_goup = remaining_tsks_godown = vFact->ntsks[0];

		srand (my_id);
		// Loop if there is some task to be processed
		while (remaining_tsks_goup > 0) {
			// Repeat the reception of messages, while there isn't waiting task
			do {	
				if (EmptyQueue(*vtsk_queue_goup)) {
					// If the queue is empty, wait until a message is waiting to be received
					flag = 1; MPI_Probe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Father, comm, &st);
				} else {
					// If the queue isn't empty, test if a message is waiting to be received
					MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Father, comm, &flag, &st);
				}
				// Receive all the pending messages
				while (flag) {
					src = st.MPI_SOURCE;
					// Receive the task and the related sizes from other processor
					MPI_Recv (vint, 3, MPI_INT, src, Tag_Send_Dims_Vector_To_Father, comm, &st);
					task = vint[0]; 
					// Define the vectors to be used in the resolution
					if (vFact[task].vecL == NULL) {
						// Create the vectors if it is necessary
						dim = vFact[task].dimX = vFact[task].dimF = vint[2];
						CreateDoubles (&vecL, 4*dim); bufL = vecL + dim; solL = bufL + dim; auxL = solL + dim;
						//Repair the structure
						vFact[task].vecL = vecL; vFact[task].bufL = bufL; 
						vFact[task].solL = solL; vFact[task].auxL = auxL; 
					} else {
						vecL = vFact[task].vecL; bufL = vFact[task].bufL; 
						solL = vFact[task].solL; auxL = vFact[task].auxL;
					}
					// Receive the vector 
					MPI_Recv (bufL+(vFact[task].dimX-vFact[task].dimF), vFact[task].dimF, MPI_DOUBLE, src,
										Tag_Send_Data_Vector_To_Father, comm, &st);
					// Verify if the father of the task is prepared to be processed
					IlupackEnQueueUp (vFact+task, vtsk_queue_goup, task);
					// Test if another message is waiting to be received
					MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Father, comm, &flag, &st);
				}
				// Try to get an actived and non-processed node of the tree
				task = DeQueue (vtsk_queue_goup, my_id);
			} while (task == -1);

			// Calculate the properties of the node
			chld = chldtab[task]; 
			// Obtain the vector on apply the preconditioner
			if (chld != -1) { // Only if the new node is not a leaf
				BuildNodeFromChildrenV (vFact, task, optP);
			}

			if (treetab[task] != -1) { // the node is not the root of the elimination tree
				// Simulation of the solution of the multilevel ILU preconditioner
				IlupackResolutionUp (vFact+task);
				// If the father is in the same node, no communication have to be done
				if (ownrtab[treetab[task]] == my_id) {
					// Verify if the father of the task is prepared to be processed
					vint[0] = task; vint[1] = vFact[task].dimX; vint[2] = vFact[task].dimF;
					IlupackEnQueueUp (vFact+task, vtsk_queue_goup, task);
				} else {
					// In other case, we have to send the data to the father
					vint[0] = task; vint[1] = vFact[task].dimX; vint[2] = vFact[task].dimF;
					smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
					MPI_Isend (vint, 3, MPI_INT, ownrtab[treetab[task]], Tag_Send_Dims_Vector_To_Father, 
											comm, &(smpnode->req));
					PushList (&lst, (void *) smpnode);
					smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
					MPI_Isend (vFact[task].bufL+(vFact[task].dimX-vFact[task].dimF), vFact[task].dimF, MPI_DOUBLE, 
											ownrtab[treetab[task]], Tag_Send_Data_Vector_To_Father, comm, &(smpnode->req));
					PushList (&lst, (void *) smpnode);
				}
			} else { // The node is the root of the elimination tree
				// Simulation of the solution of the multilevel ILU preconditioner
				IlupackResolutionUpDown (vFact+task);

				// Mark the root as processed
				IlupackEnQueueUp (vFact+task, vtsk_queue_goup, task);
				// Include the children in the godown queue
				chld = chldtab[task]; vint[0] = task; 
				if (chld != -1) {
					// The solution is copy to all the children
					BuildNodeToChildrenV (vFact, task, optP);
					// The solution is sent to all the children
					do {
						if (my_id == ownrtab[chld]) {
							// If the chld is in the same node, no communication have to be done
							EnQueue(vtsk_queue_godown, chld); 
						} else {
							smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
							MPI_Isend (&chld, 1, MPI_INT, ownrtab[chld], Tag_Send_Dims_Vector_To_Children, comm, 
													&(smpnode->req));
							PushList (&lst, (void *) smpnode);
							smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
							MPI_Isend (vFact[chld].bufL+(vFact[chld].dimX-vFact[chld].dimF), vFact[chld].dimF, MPI_DOUBLE,
													ownrtab[chld], Tag_Send_Data_Vector_To_Children, comm, &(smpnode->req));
							PushList (&lst, (void *) smpnode);
						}
						chld = brthtab[chld];
					} while (chld != -1);
				}
				// Adjust the number of processed nodes in the godown queue
       	remaining_tsks_godown--;
			}
			// Adjust the number of processed nodes in the goup queue
      	remaining_tsks_goup--;
		}
		// Remove the messages which has been arrived 
		CleanList (&lst, TestSimple);
		// Clear goup Queues if energy-shaving techniques are used 
		ClearQueue (vtsk_queue_goup);

		// Loop if there is some task to be processed
		while (remaining_tsks_godown > 0) {
			// Repeat the reception of messages, while there isn't waiting tasks
			do {	
				// Test if there is some message from some processor where a node father is
				MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Children, comm, &flag, &st);
				while (flag) {
					src = st.MPI_SOURCE;
					// Receive the task 
					MPI_Recv (&task, 1, MPI_INT, src, Tag_Send_Dims_Vector_To_Children, comm, &st);
					// Receive the vector
					MPI_Recv (vFact[task].bufL+(vFact[task].dimX-vFact[task].dimF), vFact[task].dimF, 
										MPI_DOUBLE, src, Tag_Send_Data_Vector_To_Children, comm, &st);
					// Add the received task in the queue
					EnQueue(vtsk_queue_godown, task); 
					// Test if there is some waiting task
					MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Children, comm, &flag, &st);
				}
				task = DeQueue (vtsk_queue_godown, my_id);
			} while (task == -1);

			// Simulation of the solution of the multilevel ILU preconditioner
			IlupackResolutionDown (vFact+task);

			// Include the children in the godown queue
			chld = chldtab[task];
			if (chld != -1) { // If the node is not a leaf
				// The solution is copy to all the children
				BuildNodeToChildrenV (vFact, task, optP);
				// The solution is sent to all the children
				do {
					if (my_id == ownrtab[chld]) {
						// If the chld is in the same node, no communication have to be done
						EnQueue(vtsk_queue_godown, chld); 
					} else {
						smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
						MPI_Isend (&chld, 1, MPI_INT, ownrtab[chld], Tag_Send_Dims_Vector_To_Children, comm, 
												&(smpnode->req));
						PushList (&lst, (void *) smpnode);
						smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
						MPI_Isend (vFact[chld].bufL+(vFact[chld].dimX-vFact[chld].dimF), vFact[chld].dimF, MPI_DOUBLE, 
												ownrtab[chld], Tag_Send_Data_Vector_To_Children, comm, &(smpnode->req));
						PushList (&lst, (void *) smpnode);
					}
					chld = brthtab[chld];
				} while (chld != -1);
			}
			// Adjust the number of processed nodes in the godown queue
   		remaining_tsks_godown--;
		}
		// Wait until all messages have been arrived 
		while (!emptyList (lst)) {
			CleanList (&lst, TestSimple);
		}
		// Clear goup Queues if energy-shaving techniques are used 
		ClearQueue (vtsk_queue_godown);
//	}
}

// The routine transforms to consistent, the distributed vector auxL.
void TransformVectorMPI (ptr_IlupackFactor vFact, int root, MPI_Comm comm, int optP) { 
	// Definition of the global vectors and variables
	int dimL = vFact->dimL, dimT = vFact->dimT;
	int *perm = vFact->perm, *treetab = vFact->treetab; 
	int *sizetab = vFact->auxltab;
	int *chldtab = sizetab + dimT, *nmchtab = chldtab + dimT, *brthtab = nmchtab + dimT; 
	int *wgthtab = brthtab + dimT, *marktab = wgthtab + dimT, *ownrtab = marktab + dimT;
	// Definition of the local vectors and variables
	int i, remaining_tsks = 0, remaining_tsks_goup, remaining_tsks_godown;
	double *vecL = NULL, *bufL = NULL, *solL = NULL, *auxL = NULL;
  int my_id, task = -1, chld = -1, flag, src, dim, vint[6];
	TaskQueue *vtsk_queue_goup, *vtsk_queue_godown;
	ptr_SimpleNode smpnode;
	List lst = {NULL, NULL};        // List of nonfinalized MPI_Isend 
	MPI_Status st;

	// Definition of the variable my_id
  MPI_Comm_rank(comm, &my_id);

//	if (my_id != root) {
		// Initialize the marktab vector
		InitInts (marktab, dimT, 0, 0);
		// Use the inicialize queue
		vtsk_queue_goup = vFact->goup; vtsk_queue_godown = vFact->godown;
		// Use the initial value for each queue
		ResetQueue(vtsk_queue_goup, -1);
		// Initialize the number of processes to be processed
		remaining_tsks_goup = remaining_tsks_godown = vFact->ntsks[0];
	
		srand (my_id);
		// Loop if there is some task to be processed
		while (remaining_tsks_goup > 0) {
			// Repeat the reception of messages, while there isn't waiting task
			do {	
				if (EmptyQueue(*vtsk_queue_goup)) {
					// If the queue is empty, wait until a message is waiting to be received
					flag = 1; MPI_Probe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Father, comm, &st);
				} else {
					// If the queue isn't empty, test if a message is waiting to be received
					MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Father, comm, &flag, &st);
				}
				// Receive all the pending messages
				while (flag) {
					src = st.MPI_SOURCE;
					// Receive the task and the related sizes from other processor
					MPI_Recv (vint, 3, MPI_INT, src, Tag_Send_Dims_Vector_To_Father, comm, &st);
					task = vint[0]; 
					// Define the vectors to be used in the resolution
					if (vFact[task].vecL == NULL) {
						// Create the vectors if it is necessary
						dim = vFact[task].dimX = vFact[task].dimF = vint[2];
						CreateDoubles (&vecL, 4*dim); bufL = vecL + dim; solL = bufL + dim; auxL = solL + dim;
						//Repair the structure
						vFact[task].vecL = vecL; vFact[task].bufL = bufL; 
						vFact[task].solL = solL; vFact[task].auxL = auxL; 
					} else {
						vecL = vFact[task].vecL; bufL = vFact[task].bufL; 
						solL = vFact[task].solL; auxL = vFact[task].auxL;
					}
					vFact[task].dimM = vint[2];
					// Receive the vector 
					MPI_Recv (auxL, vint[2], MPI_DOUBLE, src, Tag_Send_Data_Vector_To_Father, comm, &st);
					// Verify if the father of the task is prepared to be processed
					IlupackEnQueueUp (vFact+task, vtsk_queue_goup, task);
					// Test if another message is waiting to be received
					MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Father, comm, &flag, &st);
				}
				// Try to get an actived and non-processed node of the tree
				task = DeQueue (vtsk_queue_goup, my_id);
			} while (task == -1);

			// Calculate the properties of the node
			chld = chldtab[task]; 
			// Obtain the vector 
			if (chld != -1) { // Only if the new node is not a leaf
				CopyNodeFromChildrenV (vFact, task);
			}

			if (treetab[task] != -1) { // the node is not the root of the elimination tree
				// If the father is in the same node, no communication have to be done
				if (ownrtab[treetab[task]] == my_id) {
					// Verify if the father of the task is prepared to be processed
					vint[0] = task; vint[1] = vFact[task].dimX; vint[2] = vFact[task].dimF;
					IlupackEnQueueUp (vFact+task, vtsk_queue_goup, task);
				} else {
					// In other case, we have to send the data to the father
					vint[0] = task; vint[1] = vFact[task].dimX; vint[2] = vFact[task].dimF;
					vint[0] = task; vint[1] = vFact[task].dimM; vint[2] = vFact[task].dimM-sizetab[task];
					smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
					MPI_Isend (vint, 3, MPI_INT, ownrtab[treetab[task]], Tag_Send_Dims_Vector_To_Father, 
											comm, &(smpnode->req));
					PushList (&lst, (void *) smpnode);
					smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
					MPI_Isend (vFact[task].auxL+sizetab[task], (vFact[task].dimM-sizetab[task]), MPI_DOUBLE, 
										ownrtab[treetab[task]], Tag_Send_Data_Vector_To_Father, comm, &(smpnode->req));
					PushList (&lst, (void *) smpnode);
				}
			} else { // The node is the root of the elimination tree
				// Mark the root as processed
				IlupackEnQueueUp (vFact+task, vtsk_queue_goup, task);

				// Include the children in the godown queue
				chld = chldtab[task]; vint[0] = task; 
				if (chld != -1) {
					// The solution is sent to all the children
					do {
						if (my_id == ownrtab[chld]) {
							// If the chld is in the same node, no communication have to be done
							CopyNodeFromFatherV (vFact, chld);
							EnQueue(vtsk_queue_godown, chld); 
						} else {
							smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
							MPI_Isend (&chld, 1, MPI_INT, ownrtab[chld], Tag_Send_Dims_Vector_To_Children, 
													comm, &(smpnode->req));
							PushList (&lst, (void *) smpnode);
							smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
							MPI_Isend (vFact[task].auxL, vFact[task].dimM, MPI_DOUBLE, ownrtab[chld], 
													Tag_Send_Data_Vector_To_Children, comm, &(smpnode->req));
							PushList (&lst, (void *) smpnode);
						}
						chld = brthtab[chld];
					} while (chld != -1);
				}
				// Adjust the number of processed nodes in the godown queue
       	remaining_tsks_godown--;
			}
			// Adjust the number of processed nodes in the goup queue
     	remaining_tsks_goup--;
		}
		// Remove the messages which has been arrived 
		CleanList (&lst, TestSimple);
		// Clear goup Queues if energy-shaving techniques are used 
		ClearQueue (vtsk_queue_goup);

		// Loop if there is some task to be processed
		while (remaining_tsks_godown > 0) {
			// Repeat the reception of messages, while there isn't waiting tasks
			do {	
				// Test if there is some message from some processor where a node father is
				MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Children, comm, &flag, &st);
				while (flag) {
					src = st.MPI_SOURCE;
					// Receive the task 
					MPI_Recv (&task, 1, MPI_INT, src, Tag_Send_Dims_Vector_To_Children, comm, &st);
					// Receive the vector
					MPI_Recv (vFact[task].auxL+sizetab[task], (vFact[task].dimM-sizetab[task]), MPI_DOUBLE, 
										src, Tag_Send_Data_Vector_To_Children, comm, &st);
					// Add the received task in the queue
					EnQueue(vtsk_queue_godown, task); 
					// Test if there is some waiting task
					MPI_Iprobe (MPI_ANY_SOURCE, Tag_Send_Dims_Vector_To_Children, comm, &flag, &st);
				}
				task = DeQueue (vtsk_queue_godown, my_id);
			} while (task == -1);

			// Include the children in the godown queue
			chld = chldtab[task];
			if (chld != -1) { // If the node is not a leaf
				// The solution is sent to all the children
				do {
					vint[0] = chld; vint[1] = vFact[chld].dimX; vint[2] = vFact[chld].dimF;
					if (my_id == ownrtab[chld]) {
						// If the chld is in the same node, no communication have to be done
						CopyNodeFromFatherV (vFact, chld);
						EnQueue(vtsk_queue_godown, chld); 
					} else {
						smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
						MPI_Isend (&chld, 1, MPI_INT, ownrtab[chld], Tag_Send_Dims_Vector_To_Children, 
												comm, &(smpnode->req));
						PushList (&lst, (void *) smpnode);
						smpnode = (ptr_SimpleNode) malloc (sizeof(SimpleNode));
						MPI_Isend (vFact[task].auxL, vFact[task].dimM, MPI_DOUBLE, ownrtab[chld], 
												Tag_Send_Data_Vector_To_Children, comm, &(smpnode->req));
						PushList (&lst, (void *) smpnode);
					}
					chld = brthtab[chld];
				} while (chld != -1);
			}
			// Adjust the number of processed nodes in the godown queue
     	remaining_tsks_godown--;
		}
		// Wait until all messages have been arrived 
		while (!emptyList (lst)) {
			CleanList (&lst, TestSimple);
		}
		// Clear goup Queues if energy-shaving techniques are used 
		ClearQueue (vtsk_queue_godown);
//	}
}

// The routine simulates nItr steps of the PCG on vector rhs, obtaining the vector x,
// by using the preconditioner included in vFact. The owner of both vectors is the processor root.
// The parameter optP defines the way in which the accumulation is made.
void IlupackSolverMPI (SparseMatrix spr, double *rhs, double *sol, ptr_IlupackFactor vFact, 
												int nItr, Ilpck_Comm ilpkcomms, int optP) { 
	int i, j = 1;
  int my_id, root; 
	MPI_Comm comm, comm2;
	double te1, te2, tu1, tu2, tt1, tt2, tt;
	double scal_gbl;

	if (ilpkcomms.color) {
		// Definition of the variable my_id
		comm = ilpkcomms.comm_prec; root = ilpkcomms.root_p;
  	MPI_Comm_rank(comm, &my_id);

		MPI_Barrier (comm); reloj (&te1, &tu1); tt1 = MPI_Wtime ();
//		// Initialization of the Task Queues
//		if (vFact->ntsks == NULL) {
// 			InitTaskQueuesResolutionMPI (vFact, root, comm);
//		}
		// Initialization of the task queues and distribution of the solution
		InitResolutionMPI (vFact, rhs, root, comm);                           // vecL(NoCons)
	
//		if (my_id != root) {
		if (ilpkcomms.color2) {
			comm2 = ilpkcomms.comm_solver; 
			// Loop related with the iteration of the PCG
			for (i=0; i<nItr; i++) {
				// Simulation of the Matrix Vector Product
				if ((j == 0) && (i > 0)) SimulateProductMPI (vFact, root, comm2); // solL(Cons)   --> vecL(NoCons)
				// Application of the preconditioner
				IlupackResolutionMPI (vFact, root, comm2, optP);                  // vecL(NoCons) --> solL(Cons)
				// Simulation of the Matrix Vector Product
				if (j == 1) SimulateProductMPI (vFact, root, comm2);              // solL(Cons)   --> vecL(NoCons)
				// Global Reduction
				SimulateVectorOperationsMPI (vFact, root, comm2, 3, &scal_gbl);   // Local reduction
			}
			if (j == 1) {
				// Transform from NonConsistent to Consistent
				SimulateVectorOperationsMPI (vFact, root, comm2, 1, &scal_gbl);   // vecL(NoCons) --> auxL(NoCons)
				TransformVectorMPI (vFact, root, comm2, optP);                    // auxL(NoCons) --> auxL(Cons)
				SimulateVectorOperationsMPI (vFact, root, comm2, 2, &scal_gbl);   // aux(Cons)    --> solL(Cons)
			}
		}
	
		// Collection of the solution
		CloseResolutionMPI (vFact, sol, root, comm);
	
		reloj (&te2, &tu2); tt2 = MPI_Wtime ();
//		if (my_id == root) 
			printf ("Resolution Time (%2.2d) = %20.15e = (%20.15e,%20.15e)\n", my_id, tt2-tt1, te2-te1, tu2-tu1);
	}
}
*/

// The routine computes the multilevel ILU factorization of the sparse matrix spr.
// The initial permutation is included in the files whose names appear as parameters.
int ParmetisComputationMPI (SparseMatrix spr, int index, Ilpck_Comm ilpkcomms) {
	// Definition of the global vectors and variables
	int i, dim = spr.dim1, dimL; 
	int *perm = NULL, *rangtab = NULL, *treetab = NULL, dimT;
	int *sizetab = NULL, *chldtab = NULL, *nmchtab = NULL, *brthtab = NULL, *wgthtab = NULL;
	int remaining_tsks, *marktab = NULL, *ownrtab = NULL, *hgthtab = NULL;
	ptr_IlupackFactor vFact;
	// Definition of the local vectors and variables
	matDoubles tGlb;
	matInts mTab;
	int task = -1, chld = -1, vint[6];
  int my_id, numprocs, root; 
	int indexM = 1, ierr = 0; // JOSE
	MPI_Comm comm;
	double te1, te2, tu1, tu2, tt1, tt2, tt;
	PacketNode pcknode;
	MPI_Status st;

	// Definition of the variables numprocs and my_id
	comm = ilpkcomms.comm_metis; root = ilpkcomms.root_m;
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	dimT = 2*numprocs;
//	if (my_id == root) CreateInts (&perm, dim); 
	if (my_id == root) CreateInts (&perm, 2*dim); 
//	CreateInts (&rangtab, dimT); CreateInts (&treetab, dimT);
//	InitInts (rangtab, dimT, -1, 0); InitInts (treetab, dimT, -1, 0); 
//	CreateMatrixInts (&mTab, SIZE_MAT_TAB, dimT); InitInts (mTab[0], dimT*SIZE_MAT_TAB, 0, 0);
	CreateMatrixInts (&mTab, 2*SIZE_MAT_TAB, 2*dimT); InitInts (mTab[0], dimT*SIZE_MAT_TAB, 0, 0);
//	CreateMatrixDoubles (&tGlb, 2, SIZE_TIM_GLB); InitDoubles (tGlb[0], 2*SIZE_TIM_GLB, 0.0, 0.0);
	CreateMatrixDoubles (&tGlb, 2*2, 2*SIZE_TIM_GLB); InitDoubles (tGlb[0], 2*SIZE_TIM_GLB, 0.0, 0.0);

	rangtab = mTab[RANG]; InitInts (rangtab, dimT, -4, 0);
	treetab = mTab[TREE]; InitInts (treetab, dimT, -6, 0);
	MPI_Barrier (comm); reloj (&te1, &tu1);
	tt = ComputeMetisMPI (spr, index, root, comm, numprocs, rangtab, treetab, perm);
//	tt = ComputeMetisMPI (spr, index, root, MPI_COMM_WORLD, numprocs, rangtab, treetab, perm);
	reloj (&te2, &tu2);
	printf ("Ordering Time (%2d) = %20.15e = (%20.15e,%20.15e)\n", my_id, tt, te2-te1, tu2-tu1);
	if (my_id == root) {
		printf ("%s --> ", "rangtab"); PrintFInts (rangtab, dimT, 6, 0);
		printf ("%s --> ", "treetab"); PrintFInts (treetab, dimT, 6, 0);
	}

//	RemoveMatrixDoubles (&tGlb); RemoveMatrixInts (&mTab);
//	if (my_id == root) RemoveInts (&perm); 

	return 0;
}

/*********************************************************************************/

int main (int argc, char **argv) {
	int i, root = 0; 
	double *rhs, *sol, scale;
	SparseMatrix spr;
	DILUPACKparam paramG;
	ptr_IlupackFactor vFact = NULL;

  int my_id, numprocs;
	int rngfl = atoi(argv[2]);
	int nItr  = atoi(argv[3]);
	int nprcs = atoi(argv[4]);

	Ilpck_Comm ilpkcomms;

	int ierr = 0, index = 0, nleaves = 0;

  MPI_Init (&argc, &argv);

	// Definition of the variables numprocs and my_id
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs); MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	root = numprocs-1;
	root = 0;
	printf ("numprocs = %d , my_id = %d , root = %d\n", numprocs, my_id, root);
	// Definition of the Ilupack communicators
	CreateIlupackCommunicator (&ilpkcomms, root, nprcs);

	if (my_id == root) {
		// Read the original sparse matrix
		CreateSparseMatrixHB2 (argv[1], &spr, 1); 
		printf ("dim = %d\n", spr.dim1);
	}

	ParmetisComputationMPI (spr, index, ilpkcomms);

	printf ("numprocs = %d , my_id = %d , root = %d\n", numprocs, my_id, root);
	// Free the sparse matrix 
	if (my_id == root) {
		;
		RemoveSparseMatrix (&spr);
	}
	printf ("numprocs = %d , my_id = %d , root = %d\n", numprocs, my_id, root);
	// Free the communicators
	RemoveIlupackCommunicator (&ilpkcomms); 
//	MPI_Barrier (MPI_COMM_WORLD);
	printf ("numprocs = %d , my_id = %d , root = %d\n", numprocs, my_id, root);

	MPI_Finalize ();

	return 0;
}

