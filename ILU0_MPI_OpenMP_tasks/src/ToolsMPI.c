#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ScalarVectors.h>
#include <reloj.h>
#include <InputOutput.h>
#include "EliminationTree.h"
#include "ToolsMPI.h"
#include "Lists.h"

// #define PRINT_SEND_RESOTRNF_VECTORS 1

/*********************************************************************************/

void Sinchonization (MPI_Comm Synch_Comm, char *message) {
	int my_id, i ; 

	MPI_Comm_rank(Synch_Comm, &my_id); 
	MPI_Barrier (Synch_Comm);
	printf ("(%d) %s\n", my_id, message);
	if (my_id == 0) printf ("Waiting ... \n");
	if (my_id == 0) scanf ("%d", &i);
	if (my_id == 0) printf (" ... done\n");
	MPI_Barrier (Synch_Comm);
}

/*********************************************************************************/

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
int TestSimple (void *data) {
	int flag = 0;
	ptr_SimpleNode smpnode = (ptr_SimpleNode) data;
	
	// Verify if the communication has finalized
	MPI_Test (&(smpnode->req), &flag, &(smpnode->sta));
	if (flag) {
		// Remove the data included in the simple node
		MPI_Wait (&(smpnode->req), &(smpnode->sta));
		free (smpnode);
	}

	// Returns the result
	return flag;
}

// Detect the lost messages whose destination is one process
// into the processes of communicator Err_Comm
void DetectErrorsMPI (MPI_Comm Err_Comm) {
	int my_id, flag= 0;
	MPI_Status st;

	// Definition of the variable my_id
	MPI_Comm_rank(Err_Comm, &my_id); 
	// Test if some message exists
	MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, Err_Comm, &flag, &st);
	if (flag) {
		printf ("%d --> (%d,%d)\n", my_id, st.MPI_SOURCE, st.MPI_TAG);
	}
}

// Prepare the structures required to send/receive a ILU0Factor structure
// * Fact refers to the ILU0Factor from where the data is obtained
// * L_F defines if the matrix sprM or the matrix sprF is used
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * dimD is the size of the diagonal vectors
// * nlev is the number of level reated to the Fact
// * pcknode, where the resulted packet appears
void MakeFactorStructPacket (ILU0Factor Fact, int L_F, int size, int weight, int dimD, 
															int nlev, ptr_PacketNode pcknode) {
	int k;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
	SparseMatrix *ptr_spr = (L_F)? &(Fact.sprM): &(Fact.sprF);
		
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) ptr_spr->vptr;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size+1 ; dspl[0] = (MPI_Aint) ptr_spr->vptr ;
	type[1] = MPI_INT   ; lblq[1] = weight ; dspl[1] = (MPI_Aint) ptr_spr->vpos ;
	type[2] = MPI_DOUBLE; lblq[2] = weight ; dspl[2] = (MPI_Aint) ptr_spr->vval ;
	type[3] = MPI_DOUBLE; lblq[3] = 3*dimD ; dspl[3] = (MPI_Aint) Fact.mDia[0]  ;
	type[4] = MPI_INT   ; lblq[4] = nlev+1 ; dspl[4] = (MPI_Aint) Fact.headL    ;
	type[5] = MPI_DOUBLE; lblq[5] = nlev   ; dspl[5] = (MPI_Aint) Fact.divs     ;
	type[6] = MPI_UB    ; lblq[6] = 1      ; dspl[6] = (MPI_Aint) Fact.divs+nlev;
	for (k=6; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (7, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Prepare the structures required to send/receive a ILU0Factor structure
// * Fact refers to the ILU0Factor from where the data is obtained
// * L_F defines if the matrix sprM or the matrix sprF is used
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * nlev is the number of level reated to the Fact
// * pcknode, where the resulted packet appears
void MakeFactorStructPacket2 (ILU0Factor Fact, int L_F, int size, int weight,
																int nlev, ptr_PacketNode pcknode) {
	int k;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
	SparseMatrix *ptr_spr = (L_F)? &(Fact.sprM): &(Fact.sprF);
		
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) ptr_spr->vptr;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size+1 ; dspl[0] = (MPI_Aint) ptr_spr->vptr   ;
	type[1] = MPI_INT   ; lblq[1] = weight ; dspl[1] = (MPI_Aint) ptr_spr->vpos   ;
	type[2] = MPI_DOUBLE; lblq[2] = weight ; dspl[2] = (MPI_Aint) ptr_spr->vval   ;
	type[3] = MPI_INT   ; lblq[3] = nlev+1 ; dspl[3] = (MPI_Aint) Fact.headL    ;
	type[4] = MPI_DOUBLE; lblq[4] = nlev   ; dspl[4] = (MPI_Aint) Fact.divs     ;
	type[5] = MPI_UB    ; lblq[5] = 1      ; dspl[5] = (MPI_Aint) Fact.divs+nlev;
	for (k=5; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (6, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Prepare the structures required to send/receive a SparseMatrix structure
// * spr refers to the SparseMatrix from where the data is obtained
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * pcknode, where the resulted packet appears
void MakeSprMatrixPacket (SparseMatrix spr, int size, int weight, ptr_PacketNode pcknode) {
	int k;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
		
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) spr.vptr;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size+1; dspl[0] = (MPI_Aint) spr.vptr;
	type[1] = MPI_INT   ; lblq[1] = weight; dspl[1] = (MPI_Aint) spr.vpos;
	type[2] = MPI_DOUBLE; lblq[2] = weight; dspl[2] = (MPI_Aint) spr.vval;
	type[3] = MPI_UB; lblq[3] = 1; dspl[3] = (MPI_Aint) spr.vptr+size+1;
	for (k=3; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (4, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Prepare the structures required to send some rows of a sparse matrix 
// and its corresponding permutation value
// * spr refers to the SparseMatrix from where the data is obtained
// * src_prm includes the permutation vector
// * size is the number of used elements in src_prm
// * pcknode, where the resulted packet appears
// The routine returns the number of nonzero elements included in the packet
int MakePermSprMatrixSendPacket (SparseMatrix spr, int *src_prm, int size, 
																	ptr_PacketNode pcknode) {
	int *lblq = pcknode->lblq, *vlen = pcknode->vlen;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
	int k, row, weight = 0;

	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) src_prm;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT; lblq[0] = size; dspl[0] = (MPI_Aint) src_prm;
	type[1] = MPI_INT; lblq[1] = size; dspl[1] = (MPI_Aint) vlen;
	for (k=0; k<size; k++) { 
		// For each row of spr, some elements are defined
		row = *(src_prm+k);
		vlen[k] = spr.vptr[row+1] - spr.vptr[row];
		weight += vlen[k];
		type[k+2     ] = MPI_INT   ; lblq[k+2     ] = vlen[k]; 
		dspl[k+2     ] = (MPI_Aint) (spr.vpos+spr.vptr[row]);	
		type[k+2+size] = MPI_DOUBLE; lblq[k+2+size] = vlen[k]; 
		dspl[k+2+size] = (MPI_Aint) (spr.vval+spr.vptr[row]);
	} 
	type[2+2*size] = MPI_UB; lblq[2+2*size] = 1; 
	dspl[2+2*size] = (MPI_Aint) (vlen+size);	
	for (k=2+2*size; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (3+2*size, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));

	// Returns the number of nonzeros
	return weight;
}

// Prepare the structures required to receive some rows of a sparse matrix 
// and its corresponding permutation value
// * spr_aux refers to the SparseMatrix where the data will be stored
// * dst_prm includes the permutation vector
// * size is the number of used elements in dst_prm
// * pcknode, where the resulted packet appears
void MakePermSprMatrixRecvPacket (SparseMatrix spr_aux, int *dst_prm, int size, int weight, 
																		ptr_PacketNode pcknode) {
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
	int k;
		
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) dst_prm;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size  ; dspl[0] = (MPI_Aint) dst_prm;
	type[1] = MPI_INT   ; lblq[1] = size  ; dspl[1] = (MPI_Aint) spr_aux.vptr;
	type[2] = MPI_INT   ; lblq[2] = weight; dspl[2] = (MPI_Aint) spr_aux.vpos;
	type[3] = MPI_DOUBLE; lblq[3] = weight; dspl[3] = (MPI_Aint) spr_aux.vval;
	type[4] = MPI_UB; lblq[4] = 1; dspl[4] = (MPI_Aint) (dst_prm+size);	
	for (k=4; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (5, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
int TestPacket (void *data) {
	int flag = 0;
	ptr_PacketNode pcknode = (ptr_PacketNode) data;
	
	// Verify if the communication has finalized
	MPI_Test (&(pcknode->req), &flag, &(pcknode->sta));
	if (flag) {
		// Remove the data included in the pack
		MPI_Wait (&(pcknode->req), &(pcknode->sta));
		MPI_Type_free (&(pcknode->pack));
		free (pcknode);
	}

	// Returns the result
	return flag;
}

// From the original data (spr,perm) and the elimination tree information (treetab,rangtab),
// the routine send the data related with a leaf to prc_dst, into the communicator comm.
// In lst appears the reference to the MPI_Isend operations which are not finalized. 
// The routine returns the number of nonzero elements which are sent
int SendLeafFromMatrix_old (SparseMatrix spr, int *perm, int *treetab, int *rangtab, 
													int prc_dst, int leaf, ptr_List lst, MPI_Comm comm) {
	int i, node, weight_pack, size, weight = 0;
	ptr_PacketNode pcknode;

	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are sent by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, filled and sent
			pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
			weight_pack = MakePermSprMatrixSendPacket (spr, perm+i, size, pcknode);
			MPI_Isend (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Packet_Matrix_To_Leaf, 
									comm, &(pcknode->req));
			// The MPI_Isend operation is included in lst
			PushList (lst, (void *) pcknode);
			weight += weight_pack;
		}
		node = treetab[node];
	} while (node != -1);

	// Returns the result
	return weight;
}


// From the original matrix (spr) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst, into the communicator comm.
// In lst appears the reference to the MPI_Isend which are not finalized. 
// The routine returns the number of nonzero elements which are sent
int SendLeafFromMatrix (SparseMatrix spr, int prc_dst, int leaf, 
													ptr_ILU0Factor vFact, ptr_List lst, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int *perm = vFact->perm, *treetab = vFact->mTab[TREE], *rangtab = vFact->mTab[RANG];
	int *sizetab = vFact->mTab[SIZE], *wgthtab = vFact->mTab[WGTH];
	// Definition of the local vectors and variables
	int i, node, weight_pack, size, weight = 0;
	int dimM, nnzM, nlev, dimD;
	int *permM = NULL, *ptr = NULL, *headL = NULL, *hptr = NULL;
	double *divs = NULL;
	matDoubles mDia = NULL;
	ptr_PacketNode pcknode;

	// Calculate the properties of the node, and malloc the required vectors
	ComputeSizesNodeEliminationTree (leaf, treetab, sizetab, wgthtab, &dimM, &nnzM, &nlev);
	CreateInts (&permM, dimM+nlev+1); ptr = permM; 
	headL = permM + dimM; hptr = headL; *headL = 0;
	*(hptr+1) = *hptr + sizetab[leaf]; hptr++;
	dimD = dimM - sizetab[leaf]; 
	CreateMatrixDoubles (&mDia, SIZE_MAT_DIA, dimD); CreateDoubles (&divs , nlev); 

	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are sent by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, filled and sent
			pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
			weight_pack = MakePermSprMatrixSendPacket (spr, perm+i, size, pcknode);
			MPI_Send (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Packet_Matrix_To_Leaf, comm);
			MPI_Type_free (&(pcknode->pack));
			free (pcknode);
			weight += weight_pack;
		}
		CopyInts (perm+rangtab[node], ptr, sizetab[node]); ptr += sizetab[node];
		node = treetab[node];
	} while (node != -1);

	// Repair the structure
	vFact[leaf].dimM = dimM; vFact[leaf].dimX = dimM; vFact[leaf].permM = permM; 
	vFact[leaf].mDia = mDia; vFact[leaf].divs = divs;
	vFact[leaf].nlev = nlev; vFact[leaf].headL = headL;
	
	// Returns the result
	return weight;
}

// Compute the number of nonzeros elements of a PermSprMatrixRecvPacket packet
// * prc_src is the processor from which the messages is sent
// * sizes is the number of rows to be received
// * comm is the communicator in which the messages is sent
int ComputeSprMatrixWeights (int prc_src, int sizes, MPI_Comm comm) {
	int tam, tam_int, tam_double, tam_ub;
	MPI_Status sta;

	// Definition of sizes
	MPI_Type_size(MPI_INT, &tam_int);
	MPI_Type_size(MPI_DOUBLE, &tam_double);
	MPI_Type_size(MPI_UB, &tam_ub);
	MPI_Probe (prc_src, Tag_Send_Packet_Matrix_To_Leaf, comm, &sta);
	MPI_Get_count (&sta, MPI_BYTE, &tam);

	// Return the number of nonzeros included in a packet
	return (tam - (2*sizes*tam_int + tam_ub)) / (tam_int + tam_double);
}

// From elimination tree information (treetab,rangtab), the root sends
// the data related to a leaf:
// * sprM is the sparse matrix, and permM is the corresponding permutation.
// * headL is the vector which marks the begin and the end of each level (CSR way)
// * comm is the communicator in which the messages is received
// The routine returns the number of nonzero elements which are sent
int RecvLeafFromMatrix (int *treetab, int *rangtab, int root, int leaf,
													SparseMatrix sprM, int *permM, int *headL, MPI_Comm comm) {
	int i, nlev = 0, my_id, node, weight = 0, weight_pack, size, *pi4 = permM, dim = 0;
	int *ph = headL+1;
	SparseMatrix spr_aux;
	PacketNode pcknode;
	MPI_Status sta;

	// Definition of the variable my_id
	MPI_Comm_rank(comm, &my_id); 
	// Definition of the auxiliar sparse matrix apr_aux
	spr_aux.vptr = sprM.vptr+1; spr_aux.vpos = sprM.vpos; spr_aux.vval = sprM.vval;
	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are received by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, received and freed
			weight_pack = ComputeSprMatrixWeights (root, size, comm);
			MakePermSprMatrixRecvPacket (spr_aux, pi4, size, weight_pack, &pcknode);
			MPI_Recv (pcknode.ptr, 1, pcknode.pack, root, Tag_Send_Packet_Matrix_To_Leaf, 
									comm, &sta);
			MPI_Type_free (&(pcknode.pack));
			// The pointers are adjusted to receive the next packet
			spr_aux.vptr += size; spr_aux.vpos += weight_pack; spr_aux.vval += weight_pack;
			pi4 += size; weight += weight_pack;
		}
		*ph = rangtab[node+1]-rangtab[node]; dim += *(ph++); nlev++;
		node = treetab[node];
	} while (node != -1);
	// Transform the lenghts to the CSR way
	*(sprM.vptr) = 0; TransformLengthtoHeader (sprM.vptr, dim);
	*(headL) = 0; TransformLengthtoHeader (headL, nlev);

	// Returns the result
	return weight;
}

// The root sends the data from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
void RecvLeafFromMatrixF (ptr_ILU0Factor pFact, int root, int task, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dim, dimL = pFact->dimL;
	int *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
	int *sizetab = pFact->mTab[SIZE], *nmchtab = pFact->mTab[NMCH], *wgthtab = pFact->mTab[WGTH];
	// Definition of the local vectors and variables
	int i, dimM, indexM, dimD, nnzM, nlev;
	int *permM = NULL, *ipermM = NULL, *headL = NULL;
	double *diag = NULL, *diag1 = NULL, *diag2 = NULL, *divs = NULL;
	matDoubles mDia;
	SparseMatrix sprM;
	int numprocs, my_id;

	// Definition of the variables numprocs and my_id
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	dim = rangtab[dimL];
	// Calculate the properties of the node
	ComputeSizesNodeEliminationTree (task, treetab, sizetab, wgthtab, &dimM, &nnzM, &nlev);
	// Create the local structures
	indexM = 0; 
	CreateSparseMatrix (&sprM, indexM, dimM, dimM, nnzM, 0);
	CreateInts (&permM, dimM+nlev+1); headL = permM + dimM;
	dimD = dimM - sizetab[task];
	CreateMatrixDoubles (&mDia, SIZE_MAT_DIA, dimD); 
	diag = mDia[DIAG]; diag1 = mDia[DIA1]; diag2 = mDia[DIA2];
	CreateDoubles (&divs , nlev); CreateInts (&ipermM, dim); 
	// Fill the local structures from the shared data
	RecvLeafFromMatrix (treetab, rangtab, root, task, sprM, permM, headL, comm);
	// Eliminate incorrect elements of Sparse Matrices
	for (i=0; i<dim; i++) ipermM[i] = -1;
	ComputeInvPermutation (permM, ipermM, indexM, dimM);
	PermuteColsWithNegSparseMatrix (sprM, indexM, ipermM);
	// Readjust the size of the local data
	ReallocSparseMatrix (&sprM);
	// Scale the replicate blocks, and get diagonal before and after the scale
	GetDiagonalSparseMatrixDspls (sprM, indexM, diag1, sizetab[task]);
	AdjustLeafEliminationTree (sprM, indexM, task, treetab, nmchtab, dimL, headL, divs, nlev);
	GetDiagonalSparseMatrixDspls (sprM, indexM, diag2, sizetab[task]);
	// Compute the vector diag as (diag1-diag2)
	CopyDoubles (diag1, diag, dimD); AxpyDoubles (-1.0, diag2, diag, dimD);
	
	// Free unuseful data
	RemoveInts (&ipermM);

	// Repair the structure
	pFact->sprM  = sprM ; pFact->dimM = dimM; pFact->indM = indexM;  
	pFact->permM = permM; pFact->mDia = mDia; pFact->divs = divs ;
	pFact->nlev = nlev; pFact->headL = headL; 
}

// A slave sends the factor of task to root
// * comm is the communicator in which the messages is received
void SendFactorToMaster (ptr_ILU0Factor pFact, int task, int root, MPI_Comm comm) {
	int vint[6], *sizetab = pFact->mTab[SIZE];
	PacketNode pcknode;

	// Define the sizes to be sent to the root
	vint[0] = task; vint[1] = -1; vint[2] = -1; vint[3] = -1; vint[4] = -1; vint[5] = -1;
	if (pFact->sprF.dim1 > 0) { // If a factor has been computed
		vint[1] = pFact->dimM;
		vint[2] = pFact->sprF.dim1;
		vint[3] = pFact->sprF.vptr[vint[2]];
		vint[4] = pFact->dimM-sizetab[task];
		vint[5] = pFact->nlev;
	} 
	// Send the sizes to the root
	MPI_Send (vint, 6, MPI_INT, root, Tag_Receive_Dims_Factor_From_Leaf, comm);
	if (vint[1] != -1) { // Send the factor to the root, if it exists
		MakeFactorStructPacket (*pFact, 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
		MPI_Send (pcknode.ptr, 1, pcknode.pack, root, Tag_Receive_Data_Factor_From_Leaf, comm);
		MPI_Type_free (&(pcknode.pack));
	}
}

// The root receives a factor from a slave (src)
// * vint includes the sizes of the packet 
// * comm is the communicator in which the messages is received
int RecvFactorFromSlaveFull (ptr_ILU0Factor vFact, int src, MPI_Comm comm) {
	int vint[6];
	MPI_Status st;

	// Receive the sizes from src
	MPI_Recv (vint, 6, MPI_INT, src, Tag_Receive_Dims_Factor_From_Leaf, comm, &st);

	if (vint[1] != -1 ) {
		int indexF = 0;
		PacketNode pcknode;
		ptr_ILU0Factor pFact = vFact + vint[0];

		// Creation of a sparse matrix and other structures
		CreateSparseMatrix (&(pFact->sprF), indexF, vint[2], vint[2], vint[3], 0);
		CreateMatrixDoubles (&(pFact->mDia), SIZE_MAT_DIA, vint[4]); 
		CreateInts (&(pFact->headL) , vint[5]+1); 
		CreateDoubles (&(pFact->divs) , vint[5]); 
		// Creation of the packet and reception
		MakeFactorStructPacket (*pFact, 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
		MPI_Recv (pcknode.ptr, 1, pcknode.pack, src, Tag_Receive_Data_Factor_From_Leaf, comm, &st);
		// Destruction of the packet
		MPI_Type_free (&(pcknode.pack));
		// Adjusting of the structure
		pFact->dimX = pFact->dimF = vint[2]; pFact->nlev = vint[5];
		pFact->indF = 0; 
	}
	return vint[0];
}

// The root sends a node to a slave (src)
// * comm is the communicator in which the messages is received
// * sizeTab is the size of the leading block to factorize
void SendNodeToSlave (ptr_ILU0Factor pFact, int src, int sizeTask, MPI_Comm comm) {
	int vint[5];
	PacketNode pcknode;

	vint[0] = pFact->dimM; vint[1] = pFact->sprM.dim1;
	vint[2] = pFact->sprM.vptr[vint[1]];
	vint[3] = pFact->dimM-sizeTask;
	vint[4] = pFact->nlev;
	// Send the features of the matrix and send the data
	MPI_Send (vint, 5, MPI_INT, src, Tag_Send_Dims_Matrix_To_Leaf, comm);
	MakeFactorStructPacket (*pFact, 1, vint[1], vint[2], vint[3], vint[4], &pcknode);
	MPI_Send (pcknode.ptr, 1, pcknode.pack, src, Tag_Send_Data_Matrix_To_Leaf, comm);
	MPI_Type_free (&(pcknode.pack));
}

// A slave receives a node from the root (root)
// * comm is the communicator in which the messages is received
void RecvNodeFromMaster (ptr_ILU0Factor pFact, int root, MPI_Comm comm) {
	int indexM = 0, vint[6];
	PacketNode pcknode;
	MPI_Status st;

	// First, receive the sizes useful to create the local estructures
	MPI_Recv (vint, 6, MPI_INT, root, Tag_Send_Dims_Matrix_To_Leaf, comm, &st);
	// Update the scalars, vectors and structures, which are required to the task
	// Scalars
	pFact->dimM = vint[0]; pFact->nlev = vint[4];
	// Sparse matrix
	CreateSparseMatrix (&(pFact->sprM), indexM, vint[1], vint[1], vint[2], 0);
	// Real vectors
	CreateMatrixDoubles (&(pFact->mDia), SIZE_MAT_DIA, vint[3]); 
	CreateDoubles (&(pFact->divs) , vint[4]); 
	// Integer vectors
	CreateInts (&(pFact->permM), vint[4]+1); pFact->headL = pFact->permM;
	// Creation, reception, and freeing of the sparse matrix
	MakeFactorStructPacket (*pFact, 1, vint[1], vint[2], vint[3], vint[4], &pcknode);
	MPI_Recv (pcknode.ptr, 1, pcknode.pack, root, Tag_Send_Data_Matrix_To_Leaf, comm, &st);
	MPI_Type_free (&(pcknode.pack));
	// TO REPAIR
#ifndef _ILU0_
	pFact->indM = 0; pFact->param.testvector = NULL;
#else
	pFact->indM = 0; 
#endif
}

/*********************************************************************************/
#define NEW_ADJUST_PCG 1
// The root sends the data from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
void RecvLeafFromMatrixPCG (ptr_ILU0Factor pFact, int root, int task, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dim, dimL = pFact->dimL;
	int *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
	int *brthtab = pFact->mTab[BRTH], *chldtab = pFact->mTab[CHLD];
	int *sizetab = pFact->mTab[SIZE], *wgthtab = pFact->mTab[WGTH];
	// Definition of the local vectors and variables
	int i, dimM, indexR, nnzR, nlev;
	int *permR = NULL, *ipermR = NULL, *headL = NULL;
	SparseMatrix sprR;
	int numprocs, my_id;

	// Definition of the variables numprocs and my_id
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	dim = rangtab[dimL];
	// Calculate the properties of the node
	ComputeSizesNodeEliminationTree (task, treetab, sizetab, wgthtab, &dimM, &nnzR, &nlev);
	// Create the local structures
	indexR = 0; 
	CreateSparseMatrix (&sprR, indexR, dimM, dimM, nnzR, 0);
	permR = pFact->permM; headL = pFact->headL;
	CreateInts (&ipermR, dim);
	// Fill the local structures from the shared data
	RecvLeafFromMatrix (treetab, rangtab, root, task, sprR, permR, headL, comm);
	// Eliminate incorrect elements of Sparse Matrices
	for (i=0; i<dim; i++) ipermR[i] = -1;
	ComputeInvPermutation (permR, ipermR, indexR, dimM);
	PermuteColsWithNegSparseMatrix (sprR, indexR, ipermR);
#ifdef NEW_ADJUST_PCG
	// Remove the replicate blocks, such that only a copy exists
	AdjustLeafPCGEliminationTree (sprR, indexR, task, treetab, chldtab, brthtab, dimL, headL, nlev);
#else
	// Scale the replicate blocks, getting the diagonal before and after the operation
	AdjustLeafEliminationTree (sprR, indexR, task, treetab, nmchtab, dimL, headL, divs, nlev);
#endif
	// Readjust the size of the local data
	ReallocSparseMatrix (&sprR);
	// Free unuseful data
	RemoveInts (&ipermR);

	// Repair the structure
	pFact->sprR  = sprR ; pFact->indR = indexR;  
}

/*********************************************************************************/

// Prepare the structures required to send some elements of 
// the original vector (vec) and its corresponding permutation value.
// * vec refers to the vector from where the data is obtained
// * src_prm includes the permutation vector
// * size is the number of used elements in src_prm
// * vvec is an auxiliar vector where the data is stored before to be sent
// * pcknode, where the resulted packet appears
void MakePermVectorSendPacket (double *vec, int *src_prm, int size, double *vvec, 
																ptr_PacketNode pcknode) {
	int k, n2;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;

	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) src_prm;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size; dspl[0] = (MPI_Aint) src_prm;
	type[1] = MPI_DOUBLE; lblq[1] = size; dspl[1] = (MPI_Aint) vvec   ;
	for (k=0; k<size; k++)
		{ n2 = *(src_prm+k);  vvec[k] = vec[n2]; }
	type[2] = MPI_UB    ; lblq[2] = 1   ; dspl[2] = (MPI_Aint) (src_prm+size);	
	for (k=2; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (3, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Prepare the structures required to receive some elements of a vector
// and its corresponding permutation value
// * vec_loc refers to the vector where the data is received
// * dst_prm includes the permutation vector
// * size is the number of used elements in src_prm
// * pcknode, where the resulted packet appears
void MakePermVectorRecvPacket (double *vec_loc, int *dst_prm, int size, 
																ptr_PacketNode pcknode) {
	int k;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
		
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) dst_prm;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size; dspl[0] = (MPI_Aint) dst_prm;
	type[1] = MPI_DOUBLE; lblq[1] = size; dspl[1] = (MPI_Aint) vec_loc;
	type[2] = MPI_UB    ; lblq[2] = 1   ; dspl[2] = (MPI_Aint) (dst_prm+size);	
	for (k=2; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (3, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// From the original data (vector,perm) and the elimination tree information (treetab,rangtab),
// the routine send the data related with a leaf to prc_dst,
// where comm is the communicator in which the messages is sent.
// In lst appears the reference to the MPI_Isend which are not finalized. 
void SendLeafFromVector (double *vec, int *perm, int *treetab, int *rangtab, 
													int prc_dst, int leaf, ptr_List lst, MPI_Comm comm) {
	int i, node, size;
	ptr_PacketNode pcknode;
	double vvec[MaxPacketSize];

	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are sent by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, filled and sent
			pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
			MakePermVectorSendPacket (vec, perm+i, size, vvec, pcknode);
			MPI_Isend (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Vector_To_Leaf, 
									comm, &(pcknode->req));
			// The MPI_Isend operation is included in lst
			PushList (lst, (void *) pcknode);
		}
		node = treetab[node];
	} while (node != -1);
}

// From the original vector (vec) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst,
// where comm is the communicator in which the messages is sent.
// In lst appears the reference to the MPI_Isend which are not finalized. 
void SendLeafFromVector2 (double *vec, int prc_dst, int leaf, 
													ptr_ILU0Factor vFact, ptr_List lst, MPI_Comm comm) {
	int dimM = vFact[leaf].dimM;
	int *perm = vFact->perm, *treetab = vFact->mTab[TREE], *rangtab = vFact->mTab[RANG];
	int i, node, size;
	ptr_PacketNode pcknode;
	double *ptr = NULL;
	matDoubles mVcL;

	// Calculate the properties of the node
	CreateMatrixDoubles (&mVcL, SIZE_MAT_VCL, dimM); ptr = mVcL[VECL];
	InitDoubles (*mVcL, SIZE_MAT_VCL*dimM, -2.0, 0.0);

	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are sent by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, filled and sent
			pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
			MakePermVectorSendPacket (vec, perm+i, size, ptr, pcknode);
			MPI_Isend (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Vector_To_Leaf, 
									comm, &(pcknode->req));
			PushList (lst, (void *) pcknode);
			// The MPI_Isend operation is included in lst
			ptr += size;
		}
		node = treetab[node];
	} while (node != -1);
	
	// Repair the structure
	vFact[leaf].mVcL = mVcL;

}

// From the original vector (vec) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst,
// where comm is the communicator in which the messages is sent.
// In lst appears the reference to the MPI_Isend which are not finalized. 
void SendLeafFromVector3 (double *vec, int prc_dst, int leaf, int vid,
													ptr_ILU0Factor pFact, ptr_List lst, MPI_Comm comm) {
	int dimM = pFact->dimM;
	int *perm = pFact->perm, *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
	int i, node, size;
	ptr_PacketNode pcknode;
	double *ptr = NULL;
	double vvec[MaxPacketSize];

	// Calculate the properties of the node
	int numprocs, my_id;
	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	if (pFact->mVcL == NULL) {
		CreateMatrixDoubles (&(pFact->mVcL), SIZE_MAT_VCL, dimM);
		InitDoubles (*pFact->mVcL, SIZE_MAT_VCL*dimM, -2.0, 0.0);
		pFact->dimV = dimM; pFact->tskV = leaf;
	}
	if (pFact->mPCG == NULL) {
		CreateMatrixDoubles (&(pFact->mPCG), SIZE_MAT_PCG, pFact->dimX);
		InitDoubles (*pFact->mPCG, SIZE_MAT_PCG*pFact->dimX, -1.0, 0.0);
	}
	ptr = IdentifyVectorResolution (pFact, vid);

	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are sent by blocks
		ptr = vec + rangtab[node];
		ptr = vvec;
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, filled and sent
			pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
			MakePermVectorSendPacket (vec, perm+i, size, ptr, pcknode);
			MPI_Send (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Vector_To_Leaf, comm);
			ptr += size;
			ptr = vvec;
		}
		node = treetab[node];
	} while (node != -1);
	
}


// From elimination tree information (treetab,rangtab), the root sends
// the data related to a leaf:
// * vecL is the vector, and permM is the corresponding permutation.
// * headL is the vector which marks the begin and the end of each level (CSR way)
// * comm is the communicator in which the messages is received
void RecvLeafFromVector (int *treetab, int *rangtab, int root, int leaf,
													double *vecL, int *permM, int *headL, MPI_Comm comm) {
	int i, nlev = 0, my_id, node, size, *pi4 = permM, dim = 0;
	double *pd4 = vecL;
	int *ph = headL+1;
	PacketNode pcknode;
	MPI_Status sta;

	int numprocs;
	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	// Definition of the variable my_id
	MPI_Comm_rank(comm, &my_id); 
	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are received by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, received and freed
			MakePermVectorRecvPacket (pd4, pi4, size, &pcknode);
			MPI_Recv (pcknode.ptr, 1, pcknode.pack, root, Tag_Send_Vector_To_Leaf, 
									comm, &sta);
			MPI_Type_free (&(pcknode.pack));
			// The pointers are adjusted to receive the next packet
			pi4 += size; pd4 += size;
		}
		*ph = rangtab[node+1]-rangtab[node]; dim += *(ph++); nlev++;
		node = treetab[node];
	} while (node != -1);
	// Transform the lenghts to the CSR way
	*(headL) = 0; TransformLengthtoHeader (headL, nlev);
}

// From the original vector (vec), the root sends the corresponding elements 
// from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
void RecvLeafFromVectorF (ptr_ILU0Factor pFact, int root, int task, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dimM = pFact->dimM, nlev = pFact->nlev; 
	int *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
	int *permM = pFact->permM, *headL = pFact->headL; 
	double *divs = pFact->divs;
	double *vecL = NULL;

	if (pFact->mVcL == NULL) {
		CreateMatrixDoubles (&(pFact->mVcL), SIZE_MAT_VCL, dimM); 
		InitDoubles (*pFact->mVcL, SIZE_MAT_VCL*dimM, -2.0, 0.0);
	}
	vecL = pFact->mVcL[VECL];

	// Fill the local structures from the shared data
	RecvLeafFromVector (treetab, rangtab, root, task, vecL, permM, headL, comm);
	// Scale the vector
	AdjustVector (vecL, divs, headL, nlev);
}

// From the original vector (vec), the root sends the corresponding elements 
// from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
void RecvLeafFromVectorF3 (ptr_ILU0Factor pFact, int root, int vid, int task, int adjustVect, 
														MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dimM = pFact->dimM, nlev = pFact->nlev; 
	int *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
	int *permM = pFact->permM, *headL = pFact->headL; 
	double *divs = pFact->divs;
	double *vecL = NULL;

	int numprocs, my_id;
	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
	if (pFact->mVcL == NULL) {
		CreateMatrixDoubles (&(pFact->mVcL), SIZE_MAT_VCL, dimM);
		InitDoubles (*pFact->mVcL, SIZE_MAT_VCL*dimM, -2.0, 0.0);
		pFact->dimV = dimM; pFact->tskV = task;
	}
	if (pFact->mPCG == NULL) {
		CreateMatrixDoubles (&(pFact->mPCG), SIZE_MAT_PCG, dimM);
		InitDoubles (*pFact->mPCG, SIZE_MAT_PCG*dimM, -1.0, 0.0);
	}
	vecL = IdentifyVectorResolution (pFact, vid);

	RecvLeafFromVector (treetab, rangtab, root, task, vecL, permM, headL, comm);
	// Scale the vector
	if (adjustVect) AdjustVector (vecL, divs, headL, nlev);
}

/*********************************************************************************/

void SendVectorResolutionTransform (ptr_ILU0Factor pFact, int prc_dst, int node, int UpDown,
																		int vidI1, int vidI2, int vidO2, ptr_List lst, MPI_Comm comm) {
	int dimM = pFact->dimM, dimXX = pFact->dimX, dimF = pFact->dimF;
	int tam;
	int size;
	int *sizetab = pFact->mTab[SIZE], *chldtab = pFact->mTab[CHLD];
	double *ptr1 = NULL, *ptr2 = NULL;

	// Calculate the properties of the node
	int numprocs, my_id;
	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	MPI_Send (&size, 1, MPI_INT, prc_dst, 888, comm);
	if (pFact->mVcL == NULL) {
		printf ("ERRROOORRRR(%d)\n", my_id);
	}
	tam = (UpDown)? dimF: dimXX;
	ptr1 = IdentifyVectorResolution (pFact, vidI1) + (dimXX - tam);
#ifdef PRINT_SEND_RESOTRNF_VECTORS
	{
		int *treetab = pFact->mTab[TREE];
		char filename[80];

		sprintf (filename, "BCC_%2.2d_%2.2d", treetab[node], node);
		WriteFDoubles (filename, ptr1, tam, 40, 30);
	}
#endif
	MPI_Send (ptr1, tam, MPI_DOUBLE, prc_dst, Tag_Send_Vector_Up_1, comm);

	if (pFact->mPCG == NULL) {
		printf ("ERRROOORRRR(%d)\n", my_id);
	}
	tam = (UpDown)? sizetab[node]: 0;
  if (chldtab[node] == -1) 
		ptr2 = IdentifyVectorResolution (pFact, vidI2) + tam;
	else
		ptr2 = IdentifyVectorResolution (pFact, vidO2) + tam;
#ifdef PRINT_SEND_RESOTRNF_VECTORS
	{
		int *treetab = pFact->mTab[TREE];
		char filename[80];

		sprintf (filename, "BCD_%2.2d_%2.2d", treetab[node], node);
		WriteFDoubles (filename, ptr2, dimM-tam, 40, 30);
	}
#endif
	MPI_Send (ptr2, (dimM-tam) , MPI_DOUBLE, prc_dst, Tag_Send_Vector_Up_2, comm);
}

void RecvVectorResolutionTransform (ptr_ILU0Factor pFact, int root, int node, 
																		int vidI1, int vidI2, int vidO2, ptr_List lst, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dimM = pFact->dimM, dimXX = pFact->dimX;
	int *chldtab = pFact->mTab[CHLD];
	int tam = 0;
	double *ptr1 = NULL, *ptr2 = NULL;
	MPI_Status sta;

	int numprocs, my_id;
	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	MPI_Recv (&tam, 1, MPI_INT, root, 888, comm, &sta);
	MPI_Probe (root, Tag_Send_Vector_Up_1, comm, &sta);
	MPI_Get_count (&sta, MPI_DOUBLE, &tam);
	dimXX = tam;
	if (pFact->mVcL == NULL) {
  	CreateMatrixDoubles (&(pFact->mVcL), SIZE_MAT_VCL, dimXX);
		InitDoubles (pFact->mVcL[0], SIZE_MAT_VCL*dimXX, -3.0, 0.0);
		 pFact->dimV = dimXX; pFact->tskV = node;
	}
	if (pFact->mPCG == NULL) {
		CreateMatrixDoubles (&(pFact->mPCG), SIZE_MAT_PCG, dimM);
		InitDoubles (pFact->mPCG[0], SIZE_MAT_PCG*dimM, 0.0, 0.0);
	}
	ptr1 = IdentifyVectorResolution (pFact, vidI1);
	MPI_Recv (ptr1, dimXX, MPI_DOUBLE, root, Tag_Send_Vector_Up_1, comm, &sta);
	MPI_Get_count (&sta, MPI_DOUBLE, &tam);

	MPI_Probe (root, Tag_Send_Vector_Up_2, comm, &sta);
	MPI_Get_count (&sta, MPI_DOUBLE, &tam);
	dimM = tam;
	if (pFact->mPCG == NULL) {
		CreateMatrixDoubles (&(pFact->mPCG), SIZE_MAT_PCG, dimM);
		InitDoubles (pFact->mPCG[0], SIZE_MAT_PCG*dimM, 0.0, 0.0);
	}
	if (chldtab[node] == -1)
		ptr2 = IdentifyVectorResolution (pFact, vidI2);
	else
		ptr2 = IdentifyVectorResolution (pFact, vidO2);
	MPI_Recv (ptr2, dimM, MPI_DOUBLE, root, Tag_Send_Vector_Up_2, comm, &sta);
	pFact->dimM = dimM; pFact->dimX = dimXX;

}

/*********************************************************************************/

//
void CreateILU0Communicator (ptr_Ilpck_Comm ptr_ilpk_comm, int root, int nprocs) {
	int my_id, numprocs, my_id_n, numprocs_n;

	// Definition of the variables numprocs and my_id
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs); MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

	// Creation of the communicator for Metis computation
	MPI_Comm_dup (MPI_COMM_WORLD, &(ptr_ilpk_comm->comm_metis)); 
	// Definition of the scalars 
	ptr_ilpk_comm->root_m  = root;
	if  (root < nprocs) {
		ptr_ilpk_comm->color  = (my_id <= nprocs);
		ptr_ilpk_comm->root_p = root;
	} else {
		ptr_ilpk_comm->color  = ((my_id == root) || (my_id < nprocs));
		ptr_ilpk_comm->root_p = nprocs;
	}
	// Creation of the communicator for Preconditioner Computation
	MPI_Comm_split(ptr_ilpk_comm->comm_metis, ptr_ilpk_comm->color, my_id, &(ptr_ilpk_comm->comm_prec));
  MPI_Comm_size(ptr_ilpk_comm->comm_prec, &numprocs_n); MPI_Comm_rank(ptr_ilpk_comm->comm_prec, &my_id_n);
	// Creation of the communicator for Preconditioner Application
	ptr_ilpk_comm->color2 = ((my_id_n+1) < numprocs_n);
	MPI_Comm_split(ptr_ilpk_comm->comm_prec, ptr_ilpk_comm->color2, my_id, &(ptr_ilpk_comm->comm_solver));
  MPI_Comm_size(ptr_ilpk_comm->comm_solver, &numprocs_n); MPI_Comm_rank(ptr_ilpk_comm->comm_solver, &my_id_n);
}

void RemoveILU0Communicator (ptr_Ilpck_Comm ptr_ilpk_comm) {
	MPI_Comm_free(&(ptr_ilpk_comm->comm_solver));
	MPI_Comm_free(&(ptr_ilpk_comm->comm_prec));
	MPI_Comm_free(&(ptr_ilpk_comm->comm_metis));
	ptr_ilpk_comm->color  = ptr_ilpk_comm->color2 = -1;
	ptr_ilpk_comm->root_m = ptr_ilpk_comm->root_p = -1;
}
