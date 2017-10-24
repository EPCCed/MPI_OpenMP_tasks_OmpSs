#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <parmetis.h>
#include <ScalarVectors.h>
#include "EliminationTree.h"
#include "ToolsMPI.h"
#include "Lists.h"

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

// Prepare the structures required to send/receive a IlupackFactor structure
// * Fact refers to the IlupackFactor from where the data is obtained
// * L_F defines if the matrix sprM or the matrix sprF is used
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * dimD is the size of the diagonal vectors
// * nlev is the number of level reated to the Fact
// * pcknode, where the resulted packet appears
void MakeFactorStructPacket (IlupackFactor Fact, int L_F, int size, int weight, int dimD, 
															int nlev, ptr_PacketNode pcknode) {
	int k;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
	SparseMatrix *ptr_spr = (L_F)? &(Fact.sprM): &(Fact.sprF);
	//printf("MakeFactorStructPacket -------------------------> L_F: %d, size: %d, weight: %d, dimD: %d, nlev: %d\n", L_F, size, weight, dimD, nlev);
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) ptr_spr->vptr;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size+1 ; dspl[0] = (MPI_Aint) ptr_spr->vptr ;
	type[1] = MPI_INT   ; lblq[1] = weight ; dspl[1] = (MPI_Aint) ptr_spr->vpos ;
	type[2] = MPI_DOUBLE; lblq[2] = weight ; dspl[2] = (MPI_Aint) ptr_spr->vval ;
//	type[3] = MPI_DOUBLE; lblq[3] = 3*dimD ; dspl[3] = (MPI_Aint) Fact.diag     ;
	type[3] = MPI_DOUBLE; lblq[3] = 3*dimD ; dspl[3] = (MPI_Aint) Fact.mDia[0]  ;
	type[4] = MPI_INT   ; lblq[4] = nlev+1 ; dspl[4] = (MPI_Aint) Fact.headL    ;
	type[5] = MPI_DOUBLE; lblq[5] = nlev   ; dspl[5] = (MPI_Aint) Fact.divs     ;
	type[6] = MPI_UB    ; lblq[6] = 1      ; dspl[6] = (MPI_Aint) Fact.divs+nlev;
	for (k=6; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (7, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Prepare the structures required to send/receive a IlupackFactor structure
// * Fact refers to the IlupackFactor from where the data is obtained
// * L_F defines if the matrix sprM or the matrix sprF is used
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * nlev is the number of level reated to the Fact
// * pcknode, where the resulted packet appears
void MakeFactorStructPacket2 (IlupackFactor Fact, int L_F, int size, int weight,
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
//	printf ("MakePermSprMatrixSendPacket - BEFORE MPI_Type_create_struct\n");
	MPI_Type_create_struct (3+2*size, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
//	printf ("MakePermSprMatrixSendPacket - AFTER  MPI_Type_commit\n");

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
//	printf ("MakePermSprMatrixRecvPacket - BEFORE MPI_Type_create_struct (%d,%d)\n", size, weight);
	MPI_Type_create_struct (5, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
//	printf ("MakePermSprMatrixRecvPacket - AFTER  MPI_Type_commit\n");
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

#define SENDLEAF_ISEND 1 

// From the original matrix (spr) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst, into the communicator comm.
// In lst appears the reference to the MPI_Isend which are not finalized. 
// The routine returns the number of nonzero elements which are sent
int SendLeafFromMatrix (SparseMatrix spr, int prc_dst, int leaf, 
													ptr_IlupackFactor vFact, ptr_List lst, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dimT = vFact->dimT;
//	int *perm = vFact->perm, *treetab = vFact->treetab, *rangtab = vFact->rangtab;
	int *perm = vFact->perm, *treetab = vFact->mTab[TREE], *rangtab = vFact->mTab[RANG];
//	int *sizetab = vFact->auxltab, *wgthtab = sizetab + 4*dimT;
	int *sizetab = vFact->mTab[SIZE], *wgthtab = vFact->mTab[WGTH];
	// Definition of the local vectors and variables
	int i, node, weight_pack, size, weight = 0;
	int dimM, nnzM, nlev, dimD;
	int *permM = NULL, *ptr = NULL, *headL = NULL, *hptr = NULL;
//	double *diag = NULL;
	double *divs = NULL;
	matDoubles mDia = NULL;
	ptr_PacketNode pcknode;

	// Calculate the properties of the node, and malloc the required vectors
	ComputeSizesNodeEliminationTree (leaf, treetab, sizetab, wgthtab, &dimM, &nnzM, &nlev);
	CreateInts (&permM, dimM+nlev+1); ptr = permM; 
	headL = permM + dimM; hptr = headL; *headL = 0;
//	*(hptr+1) = *hptr + sizetab[node]; hptr++;
	*(hptr+1) = *hptr + sizetab[leaf]; hptr++;
//	dimD = dimM - sizetab[leaf]; CreateDoubles (&diag, 3*dimD+nlev);
	dimD = dimM - sizetab[leaf]; 
	CreateMatrixDoubles (&mDia, SIZE_MAT_DIA, dimD); CreateDoubles (&divs , nlev); 

	// The elimination tree is walking from the leaf to the root, 
	node = leaf;
	do {
		// The rows and corresponding permutation values are sent by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
//			printf ("(%d) Sending row %d to %d\n", 0, i, prc_dst);
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, filled and sent
			pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
//			printf ("(%d) Making packet of row %d (%d) to %d\n", 0, i, size, prc_dst);
			weight_pack = MakePermSprMatrixSendPacket (spr, perm+i, size, pcknode);
//			printf ("(%d) Made packet of row %d (%d,%d) to %d\n", 0, i, size, weight_pack, prc_dst);
#ifdef SENDLEAF_ISEND
			MPI_Isend (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Packet_Matrix_To_Leaf, 
									comm, &(pcknode->req));
			// The MPI_Isend operation is included in lst
			PushList (lst, (void *) pcknode);
#else
			MPI_Send (pcknode->ptr, 1, pcknode->pack, prc_dst, Tag_Send_Packet_Matrix_To_Leaf, comm);
//			printf ("(%d) Sent row %d to %d\n", 0, i, prc_dst);
			MPI_Type_free (&(pcknode->pack));
			free (pcknode);
#endif
			weight += weight_pack;
		}
		CopyInts (perm+rangtab[node], ptr, sizetab[node]); ptr += sizetab[node];
		node = treetab[node];
	} while (node != -1);

	// Repair the structure
	vFact[leaf].dimM = dimM; vFact[leaf].permM = permM; 
	vFact[leaf].mDia = mDia; vFact[leaf].divs  = divs ;
//	vFact[leaf].diag  = diag                    ; vFact[leaf].diag1 = vFact[leaf].diag  + dimD; 
//	vFact[leaf].diag2 = vFact[leaf].diag1 + dimD; vFact[leaf].divs  = vFact[leaf].diag2 + dimD; 
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
	//printf("0 RecvLeafFromMatrix\n");
	do {
		// The rows and corresponding permutation values are received by blocks
		for (i=rangtab[node]; i<rangtab[node+1]; i+=MaxPacketSize) {
//			printf ("(%d) Receiving row %d\n", my_id, i);
			// The number of rows to sent is defined
			size = rangtab[node+1] - i;
	//		printf("1 RecvLeafFromMatrix\n");
			if (size > MaxPacketSize) size = MaxPacketSize;
			// The packet is created, received and freed
			weight_pack = ComputeSprMatrixWeights (root, size, comm);
		//	printf ("(%d) Making packet of row %d (%d,%d)\n", my_id, i, size, weight_pack);
			MakePermSprMatrixRecvPacket (spr_aux, pi4, size, weight_pack, &pcknode);
		//	printf ("(%d) Made packet of row %d (%d,%d) from %d\n", my_id, i, size, weight_pack, root);
			MPI_Recv (pcknode.ptr, 1, pcknode.pack, root, Tag_Send_Packet_Matrix_To_Leaf, 
									comm, &sta);
		//	printf ("(%d) Received row %d\n", my_id, i);
			MPI_Type_free (&(pcknode.pack));
			// The pointers are adjusted to receive the next packet
			spr_aux.vptr += size; spr_aux.vpos += weight_pack; spr_aux.vval += weight_pack;
			pi4 += size; weight += weight_pack;
		}
		*ph = rangtab[node+1]-rangtab[node]; dim += *(ph++); nlev++;
		node = treetab[node];
	} while (node != -1);
//	printf("RecvLeafFromMatrix\n");
	// Transform the lenghts to the CSR way
	*(sprM.vptr) = 0; TransformLengthtoHeader (sprM.vptr, dim);
	*(headL) = 0; TransformLengthtoHeader (headL, nlev);

	// Returns the result
	return weight;
}

// The root sends the data from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
void RecvLeafFromMatrixF (ptr_IlupackFactor pFact, int root, int task, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dim, dimL = pFact->dimL, dimT = pFact->dimT;
//	int *treetab = pFact->treetab, *rangtab = pFact->rangtab;
	int *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
//	int *sizetab = pFact->auxltab, *nmchtab = sizetab + 2*dimT, *wgthtab = sizetab + 4*dimT;
	int *sizetab = pFact->mTab[SIZE], *nmchtab = pFact->mTab[NMCH], *wgthtab = pFact->mTab[WGTH];
	// Definition of the local vectors and variables
	int i, dimM, indexM, dimD, nnzM, nlev;
	int *permM = NULL, *ipermM = NULL, *headL = NULL;
	double *diag = NULL, *diag1 = NULL, *diag2 = NULL, *divs = NULL;
	matDoubles mDia;
	SparseMatrix sprM;
	int numprocs, my_id;
//	char nameFile[80];

	// Definition of the variables numprocs and my_id
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

//	dim = rangtab[dimT-1];
	dim = rangtab[dimL];
	// Calculate the properties of the node
	ComputeSizesNodeEliminationTree (task, treetab, sizetab, wgthtab, &dimM, &nnzM, &nlev);
//	printf ("(%d) dimM = %d , nnzM = %d , nlev = %d\n", my_id, dimM, nnzM, nlev);
	// Create the local structures
	indexM = 0; // JOSE
	CreateSparseMatrix (&sprM, indexM, dimM, dimM, nnzM, 0);
	CreateInts (&permM, dimM+nlev+1); headL = permM + dimM;
	dimD = dimM - sizetab[task];
	CreateMatrixDoubles (&mDia, SIZE_MAT_DIA, dimD); 
	diag = mDia[DIAG]; diag1 = mDia[DIA1]; diag2 = mDia[DIA2];
	CreateDoubles (&divs , nlev); CreateInts (&ipermM, dim); 
	// Fill the local structures from the shared data
//	printf ("(%d) BEFORE RecvLeafFromMatrix, task = %d\n", my_id, task);
	RecvLeafFromMatrix (treetab, rangtab, root, task, sprM, permM, headL, comm);
//	printf ("(%d) AFTER  RecvLeafFromMatrix, task = %d\n", my_id, task);
//	sprintf (nameFile, "MatrixA_%2.2d.rsa", task); WriteSparseMatrixHB (nameFile, sprM, 10, 4, "ILU", 1);
	// Eliminate incorrect elements of Sparse Matrices
	for (i=0; i<dim; i++) ipermM[i] = -1;
	ComputeInvPermutation (permM, ipermM, indexM, dimM);
	PermuteColsWithNegSparseMatrix (sprM, indexM, ipermM);
//	sprintf (nameFile, "MatrixB_%2.2d.rsa", task); WriteSparseMatrixHB (nameFile, sprM, 10, 4, "ILU", 1);
	// Readjust the size of the local data
	ReallocSparseMatrix (&sprM);
	// Scale the replicate blocks, and get diagonal before and after the scale
	GetDiagonalSparseMatrixDspls (sprM, indexM, diag1, sizetab[task]);
	AdjustLeafEliminationTree (sprM, indexM, task, treetab, nmchtab, dimL, headL, divs, nlev);
//	sprintf (nameFile, "MatrixC_%2.2d.rsa", task); WriteSparseMatrixHB (nameFile, sprM, 10, 4, "ILU", 1);
	GetDiagonalSparseMatrixDspls (sprM, indexM, diag2, sizetab[task]);
	// Compute the vector diag as (diag1-diag2)
	CopyDoubles (diag1, diag, dimD); AxpyDoubles (-1.0, diag2, diag, dimD);
	
	// Free unuseful data
	RemoveInts (&ipermM);

	// Repair the structure
	pFact->sprM  = sprM ; pFact->dimM = dimM; pFact->indM = indexM;  
	pFact->permM = permM; pFact->mDia = mDia; pFact->divs = divs ;
//	pFact->diag = diag; pFact->diag1 = diag1; pFact->diag2 = diag2;
	pFact->nlev = nlev; pFact->headL = headL; 
	// pFact->divs  = divs ;
}

// A slave sends the factor of task to root
// * comm is the communicator in which the messages is received
void SendFactorToMaster (ptr_IlupackFactor pFact, int task, int root, MPI_Comm comm) {
	int indexF = 0, vint[6], *sizetab = pFact->mTab[SIZE];
	PacketNode pcknode;
	int err;
	//printf("***************DENTRO SendFactorToMaster ************ task: %d, nodo_al_que_envio: %d, pFact: %d\n", task, root, pFact);
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
  //printf("******Send en SndFactor: "); PrintInts (vint, 6);
	if (vint[1] != -1) { // Send the factor to the root, if it exists
		MakeFactorStructPacket (*pFact, 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
		err = MPI_Send (pcknode.ptr, 1, pcknode.pack, root, Tag_Receive_Data_Factor_From_Leaf, comm);
//		MPI_Barrier (MPI_COMM_WORLD); MPI_Finalize (); exit(0);
		MPI_Type_free (&(pcknode.pack));
	}
//	printf("FIN DE SEND FACTOR TO MASTER\n");
}

// The root receives a factor from a slave (src)
// * vint includes the sizes of the packet 
// * comm is the communicator in which the messages is received
void RecvFactorFromSlave (ptr_IlupackFactor pFact, int *vint, int src, MPI_Comm comm) {
	int indexF = 0;
	PacketNode pcknode;
	MPI_Status st;
  PrintInts (vint, 6);
	// Creation of a sparse matrix and other structures
	CreateSparseMatrix (&(pFact->sprF), indexF, vint[2], vint[2], vint[3], 0);
	CreateMatrixDoubles (&(pFact->mDia), SIZE_MAT_DIA, vint[4]); 
	CreateInts (&(pFact->headL), vint[5]+1); 
	CreateDoubles (&(pFact->divs) , vint[5]); 
	// Creation of the packet and reception
	MakeFactorStructPacket (*pFact, 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
	MPI_Recv (pcknode.ptr, 1, pcknode.pack, src, Tag_Receive_Data_Factor_From_Leaf, comm, &st);
	//printf("3--***************DENTRO RecvFactorFromSlave_New ************\n");
	// Destruction of the packet
	MPI_Type_free (&(pcknode.pack));
	// Adjusting of the structure
	pFact->dimF = vint[2]; pFact->nlev = vint[5];
	// TO REPAIR
	pFact->indF = 1; pFact->param.testvector = NULL;
}

// The root receives a factor from a slave (src)
// * comm is the communicator in which the messages is received
void RecvFactorFromSlave_New (ptr_IlupackFactor pFact, int src, MPI_Comm comm) {
	int indexF = 0, vint[6];
	PacketNode pcknode;
	MPI_Status st;
	int err; 
  //printf("******Entro en RcvFactor, pFact: %d\n", pFact);
	// Receive the sizes to the root
	MPI_Recv (vint, 6, MPI_INT, src, Tag_Receive_Dims_Factor_From_Leaf, comm, &st);
  //printf("******Recv en RcvFactor (1): "); PrintInts (vint, 6);
	if (vint[1] != -1) { //LO HE PUESTO YO
		// Creation of a sparse matrix and other structures
		CreateSparseMatrix (&(pFact->sprF), indexF, vint[2], vint[2], vint[3], 0);
		CreateMatrixDoubles (&(pFact->mDia), SIZE_MAT_DIA, vint[4]); 
		CreateInts (&(pFact->headL), vint[5]+1); 
		CreateDoubles (&(pFact->divs) , vint[5]); 
  	//printf("******Recv en RcvFactor (2): "); PrintInts (vint, 6);
		//printf("1--***************DENTRO RecvFactorFromSlave_New ************\n");
		// Creation of the packet and reception
		MakeFactorStructPacket (*pFact, 0, vint[2], vint[3], vint[4], vint[5], &pcknode);
		//printf("2--***************DENTRO RecvFactorFromSlave_New ************, src: %d \n", src);
		err = MPI_Recv (pcknode.ptr, 1, pcknode.pack, src, Tag_Receive_Data_Factor_From_Leaf, comm, &st);
//		MPI_Barrier (MPI_COMM_WORLD); MPI_Finalize (); exit(0);
		// Destruction of the packet
		MPI_Type_free (&(pcknode.pack));
		// Adjusting of the structure
		pFact->dimF = vint[2]; pFact->nlev = vint[5];
		// TO REPAIR
		pFact->indF = 1; pFact->param.testvector = NULL;
	}
//	printf("FIN DE RECV FACTOR FROM SLAVE\n");
}

// The root sends a node to a slave (src)
// * comm is the communicator in which the messages is received
// * sizeTab is the size of the leading block to factorize
void SendNodeToSlave (ptr_IlupackFactor pFact, int src, int sizeTask, MPI_Comm comm) {
	int vint[5];
	PacketNode pcknode;

	vint[0] = pFact->dimM; vint[1] = pFact->sprM.dim1;
	vint[2] = pFact->sprM.vptr[vint[1]];
	vint[3] = pFact->dimM-sizeTask;
	vint[4] = pFact->nlev;
	// Send the features of the matrix and send the data
	MPI_Send (vint, 5, MPI_INT, src, Tag_Send_Dims_Matrix_To_Leaf, comm);
//	MakeFactorStructPacket2 (vFact[task], 1, vint[1], vint[2], vint[4], &pcknode);
	MakeFactorStructPacket (*pFact, 1, vint[1], vint[2], vint[3], vint[4], &pcknode);
	MPI_Send (pcknode.ptr, 1, pcknode.pack, src, Tag_Send_Data_Matrix_To_Leaf, comm);
	MPI_Type_free (&(pcknode.pack));
}

// A slave receives a node from the root (root)
// * comm is the communicator in which the messages is received
void RecvNodeFromMaster (ptr_IlupackFactor pFact, int root, MPI_Comm comm) {
	//printf("---------------------------------------------------------------------INICIO RCV NODE FROM MASTER. NODO DE RECEPCION: %d", root);
	int indexM = 0, vint[6];
	PacketNode pcknode;
	MPI_Status st;
	//printf("INICIO RCV NODE FROM MASTER. NODO DE RECEPCION: %d", root);
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
	pFact->indM = 0; pFact->param.testvector = NULL;
	//printf("FIN RCV NODE FROM MASTER");
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
		{ n2 = *(src_prm+k); /* n2 = perm[n1]; */ vvec[k] = vec[n2]; } 
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
													ptr_IlupackFactor vFact, ptr_List lst, MPI_Comm comm) {
//	int dimT = vFact->dimT, dimM = vFact[leaf].dimM;
	int dimM = vFact[leaf].dimM;
//	int *perm = vFact->perm, *treetab = vFact->treetab, *rangtab = vFact->rangtab;
	int *perm = vFact->perm, *treetab = vFact->mTab[TREE], *rangtab = vFact->mTab[RANG];
//	int *sizetab = vFact->auxltab; // , *wgthtab = sizetab + 4*dimT;
	int i, node, size;
	ptr_PacketNode pcknode;
	double *vecL = NULL, *ptr = NULL;
	matDoubles mVcL;
//	char nameFile[80];

	// Calculate the properties of the node
//	CreateDoubles (&vecL, 4*dimM); ptr = vecL;
	CreateMatrixDoubles (&mVcL, SIZE_MAT_VCL, dimM); ptr = mVcL[VECL];

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
//	vFact[leaf].vecL = vecL;                    vFact[leaf].bufL = vFact[leaf].vecL + dimM; 
//	vFact[leaf].solL = vFact[leaf].bufL + dimM; vFact[leaf].auxL = vFact[leaf].solL + dimM;

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
void RecvLeafFromVectorF (ptr_IlupackFactor pFact, int root, int task, MPI_Comm comm) {
	// Definition of the global vectors and variables
	int dimM = pFact->dimM, nlev = pFact->nlev; 
//	int *treetab = pFact->treetab, *rangtab = pFact->rangtab;
	int *treetab = pFact->mTab[TREE], *rangtab = pFact->mTab[RANG];
	int *permM = pFact->permM, *headL = pFact->headL; 
	double *divs = pFact->divs;
	double *vecL = NULL;
//	double *vecL = pFact->vecL, *bufL = pFact->bufL;
//	double *solL = pFact->solL, *auxL = pFact->auxL;

	// Create the local structures
//	if (pFact->vecL == NULL) {
//		CreateDoubles (&vecL, 4*dimM); bufL = vecL + dimM; solL = bufL + dimM; auxL = solL + dimM;
//		//Repair the structure
//		pFact->vecL = vecL; pFact->bufL = bufL; pFact->solL = solL; pFact->auxL = auxL; 
//	} else
//		vecL = pFact->vecL;
	if (pFact->mVcL == NULL)
		CreateMatrixDoubles (&(pFact->mVcL), SIZE_MAT_VCL, dimM); 
	vecL = pFact->mVcL[VECL];

	// Fill the local structures from the shared data
	RecvLeafFromVector (treetab, rangtab, root, task, vecL, permM, headL, comm);
	// Scale the vector
	AdjustVector (vecL, divs, headL, nlev);
}

/*********************************************************************************/

//
void CreateIlupackCommunicator (ptr_Ilpck_Comm ptr_ilpk_comm, int root, int nprocs) {
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
	if (ptr_ilpk_comm->color)
//	printf ("%2.2d -> (PREC) %2.2d de %2.2d con root = %2.2d\n", my_id, my_id_n, numprocs_n, ptr_ilpk_comm->root_p);
	// Creation of the communicator for Preconditioner Application
	ptr_ilpk_comm->color2 = ((my_id_n+1) < numprocs_n);
//	MPI_Comm_split(ptr_ilpk_comm->comm_prec, ((my_id_n+1) == numprocs_n), my_id, &(ptr_ilpk_comm->comm_solver));
	MPI_Comm_split(ptr_ilpk_comm->comm_prec, ptr_ilpk_comm->color2, my_id, &(ptr_ilpk_comm->comm_solver));
  MPI_Comm_size(ptr_ilpk_comm->comm_solver, &numprocs_n); MPI_Comm_rank(ptr_ilpk_comm->comm_solver, &my_id_n);
//	if (ptr_ilpk_comm->color)
//	printf ("%2.2d -> (SOLV) %2.2d de %2.2d\n", my_id, my_id_n, numprocs_n);
}

//
void RemoveIlupackCommunicator (ptr_Ilpck_Comm ptr_ilpk_comm) {
	MPI_Comm_free(&(ptr_ilpk_comm->comm_solver));
	MPI_Comm_free(&(ptr_ilpk_comm->comm_prec));
	MPI_Comm_free(&(ptr_ilpk_comm->comm_metis));
	ptr_ilpk_comm->color  = ptr_ilpk_comm->color2 = -1;
	ptr_ilpk_comm->root_m = ptr_ilpk_comm->root_p = -1;
}

//
void ComputeTreeTab (int *treetab, int numprocs) {
	int numL = numprocs, numN = 2*numL-1;
	int i, j = 0, k = numL;

	while (numL > 1) {
		for (i=0; i<numL; i+=2) {
			treetab[i+j] = k;
			treetab[i+j+1] = k;
//			printf ("[%d] = %d\n", i+j  , k);
//			printf ("[%d] = %d\n", i+j+1, k);
			k++;
		}
		j += numL; numL /= 2;
	}
	treetab[j] = -1;
//	printf ("[%d] = %d\n", j, -1);
}

//
double ComputeMetisMPI (SparseMatrix spr, int index, int root, MPI_Comm comm_m,
													int nleaves, int *rangtab, int *treetab, int *perm) {
//	int i, dim = *Dim, dim_local;
	int i, dim = spr.dim1, dim_local, size_local;
	int my_id, numprocs;
	int *vptr = NULL, *vpos = NULL;
	int *vptr_loc = NULL, *vpos_loc = NULL;
	int *sizes = NULL, *weights = NULL, *dspls = NULL;
	int *iperm = NULL, *sizes2 = NULL;
	int numflag;
	int indexM; // JOSE
	int indexG = index; // JOSE
	int *order = NULL, *options = NULL;
	double t1, t2, t_parm;

	MPI_Comm_size(comm_m, &numprocs);	MPI_Comm_rank(comm_m, &my_id);

/*********** COMPUTATION OF THE GRAPH  ***************/

	if (my_id == root) {
		// Compute the direct graph from undirect graph
		CreateInts (&vptr, dim+1); CreateInts (&vpos, spr.vptr[dim]*2);
		GetGraphSparseMatrix (spr, index, vptr, vpos, indexG);
	}

/*********** DISTRIBUTION OF THE GRAPH  ***************/

	// Malloc integer vectors for ParMetis
	if (my_id == root) CreateInts (&iperm, dim); 
	CreateInts (&sizes2, 2*nleaves); 

	// Malloc integer vectors for Distribution
	CreateInts (&sizes, numprocs); CreateInts (&weights, numprocs);
	CreateInts (&dspls, numprocs+1); 

	// Compute the number of rows in each processor
	MPI_Bcast (&dim, 1, MPI_INT, root, comm_m); 
	for (i=0; i < numprocs; i++) sizes[i] = (dim / numprocs);
	for (i=0; i < (dim % numprocs); i++) sizes[i]++; 

	// Distribute the sizes of the rows
	dim_local = sizes[my_id]; 
	if (my_id == root) {
		// Obtain the vlen vector from vptr vector
		TransformHeadertoLength (vptr, dim);

//		CopyInts (sizes, dspls+1, numprocs); *dspls = 0; TransformLengthtoHeader (dspls, numprocs);
		*dspls = 0; ComputeHeaderfromLength (sizes, dspls, numprocs);

		for (i=0; i<numprocs; i++) weights[i] = AddInts (vptr+1+dspls[i], dspls[i+1] - dspls[i]);
	}
//	else {
//		CreateInts (&vptr, dim_local+1);
//	}
		CreateInts (&vptr_loc, dim_local+1);

//	MPI_Scatterv (vptr+1, sizes, dspls, MPI_INT, (my_id == root)?(vptr+1+dspls[root]):(vptr+1), 
	MPI_Scatterv (vptr+1, sizes, dspls, MPI_INT, (vptr_loc+1), 
								dim_local, MPI_INT, root, comm_m);
//								dim_local, MPI_INT, root, MPI_COMM_WORLD);

	// Compute the number of nonzeros in each processor
	if (my_id == root) {
		*dspls = 0; ComputeHeaderfromLength (weights, dspls, numprocs);
	}

	// Distribute the position of the nonzeros 
//	MPI_Scatter (weights, 1, MPI_INT, &size_local, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Scatter (weights, 1, MPI_INT, &size_local, 1, MPI_INT, root, comm_m);
//	if (my_id != 0) CreateInts (&vpos, size_local);
	CreateInts (&vpos_loc, size_local);
	// Distribute the nonzero patterns
//	MPI_Scatterv (vpos, weights, dspls, MPI_INT, (my_id == root)?MPI_IN_PLACE:vpos, 
	MPI_Scatterv (vpos, weights, dspls, MPI_INT, vpos_loc, 
								size_local, MPI_INT, root, comm_m);
//								size_local, MPI_INT, root, MPI_COMM_WORLD);

	// Obtain the vptr vector from len vector
	*vptr_loc = 0; TransformLengthtoHeader (vptr_loc, dim_local);

	CopyInts (sizes, dspls+1, numprocs); *dspls = 0; TransformLengthtoHeader (dspls, numprocs);

/*********** COMPUTATION OF THE ELIMINATION GRAPH  ***************/

	// Compute the ParMetis permutation
	numflag = 0; 
//	CreateInts (&order, dim_local); for (i=0; i<dim_local; i++) order[i] = dspls[my_id] + i;
	CreateInts (&order, (dim_local+1)*numprocs); for (i=0; i<dim_local; i++) order[i] = dspls[my_id] + i;
	CreateInts (&options, 3); options[0] = 3; options[1] = 0; options[2] = 1;
	MPI_Barrier (comm_m); t1 = MPI_Wtime ();
	ParMETIS_V3_NodeND(dspls, vptr_loc, vpos_loc, &numflag, options, order, sizes2, &comm_m);
//                                         vwgt           mtype rtype pnsps snsps ubfrc  seed dbglvl
//	ParMETIS_V32_NodeND(dspls, vptr, vpos, NULL, &numflag, NULL, NULL, NULL, NULL, NULL, NULL, NULL, order, sizes2, &comm);
	t2 = MPI_Wtime (); t_parm = t2 - t1; 

	*rangtab = 0; ComputeHeaderfromLength (sizes2, rangtab, 2*nleaves-1);

	// Distribute the number of rows of each processor from the elimination tree
	MPI_Gatherv (order, dim_local, MPI_INT, iperm, sizes, dspls, MPI_INT, root, comm_m);

	ComputeTreeTab (treetab, numprocs);

	// Free Parmetis and Graph Structures
	RemoveInts (&options); RemoveInts (&order); 

	// Compute the number of nonzeros of each node of the elimination tree
	if (my_id == root) {
		// Obtain the distribution permutation
		indexM = 0; // JOSE
		ComputeInvPermutation (iperm, perm, indexM, dim);
	}

	RemoveInts (&vpos_loc); RemoveInts (&vptr_loc);
	RemoveInts (&dspls); RemoveInts (&weights); RemoveInts (&sizes); 
	RemoveInts (&sizes2); 
	if (my_id == root) {
		RemoveInts (&iperm); 
		RemoveInts (&vpos); RemoveInts (&vptr);
	}

	return t_parm;
}

