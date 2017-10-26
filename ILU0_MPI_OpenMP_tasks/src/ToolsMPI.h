#ifndef ToolsMPI

#define ToolsMPI 1

#include "ILU0Factor.h"
#include "Lists.h"

/*********************************************************************************/

#define Tag_Demand_Matrix_From_Root       1001
#define Tag_Send_Task_To_Leaf             1002
#define Tag_Receive_Dims_Factor_From_Leaf 1003
#define Tag_End_Distribution_To_Leaf      1004
#define Tag_Send_Dims_Matrix_To_Leaf      1006
#define Tag_Send_Data_Matrix_To_Leaf      1007

#define Tag_Demand_Vector_From_Root       1011
#define Tag_Send_Dims_Vector_To_Father    1015
#define Tag_Send_Data_Vector_To_Father    1016

#define Tag_Send_Task_To_Root             1021
#define Tag_Send_Solution_To_Root         1022
#define Tag_Send_Dims_Vector_To_Children  1025
#define Tag_Send_Data_Vector_To_Children  1026

#define Tag_End_Resolution_To_Leaf        1031

#define Tag_Send_Vector_Up_1              1041
#define Tag_Send_Vector_Up_2              1042

#define Tag_Send_Packet_Matrix_To_Leaf     210
#define Tag_Receive_Data_Factor_From_Leaf  220
#define Tag_Send_Vector_To_Leaf            230

/*********************************************************************************/

// typedef struct SimpleNode {
typedef struct {
	MPI_Status sta;
	MPI_Request req;
} SimpleNode, *ptr_SimpleNode;

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
extern int TestSimple (void *data);

/*********************************************************************************/

#define MaxPacketSize                    5000

typedef struct {
	unsigned char *ptr;
	int lblq[2*MaxPacketSize+3], vlen[MaxPacketSize];
	MPI_Aint dspl[2*MaxPacketSize+3];
	MPI_Datatype type[2*MaxPacketSize+3];
	MPI_Datatype pack;
	MPI_Status sta;
	MPI_Request req;
} PacketNode, *ptr_PacketNode;

/*********************************************************************************/

extern void Sinchonization (MPI_Comm Synch_Comm, char *message);

/*********************************************************************************/

// Detect the lost messages whose destination is one process
// into the processes of communicator Err_Comm
extern void DetectErrorsMPI (MPI_Comm Err_Comm);

// Prepare the structures required to send/receive a ILU0Factor structure
// * Fact refers to the ILU0Factor from where the data is obtained
// * L_F defines if the matrix sprL or the matrix sprF is used
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * dimD is the size of the diagonal vectors
// * nlev is the number of level reated to the Fact
// * pcknode, where the result appears
extern void MakeFactorStructPacket (ILU0Factor Fact, int L_F, int size, int weight, int dimD, 
																			int nlev, ptr_PacketNode pcknode);

// Prepare the structures required to send/receive a ILU0Factor structure
// * Fact refers to the ILU0Factor from where the data is obtained
// * L_F defines if the matrix sprL or the matrix sprF is used
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * nlev is the number of level reated to the Fact
// * pcknode, where the result appears
extern void MakeFactorStructPacket2 (ILU0Factor Fact, int L_F, int size, int weight,
																			int nlev, ptr_PacketNode pcknode);

// Prepare the structures required to send/receive a SparseMatrix structure
// * spr refers to the SparseMatrix from where the data is obtained
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * pcknode, where the result appears
extern void MakeSprMatrixPacket (SparseMatrix spr, int size, int weight, ptr_PacketNode pcknode);

// Prepare the structures required to send some rows of a sparse matrix 
// and its corresponding permutation value
// * spr refers to the SparseMatrix from where the data is obtained
// * src_prm includes the permutation vector
// * size is the number of used elements in src_prm
// * pcknode, where the result appears
// The routine returns the number of nonzero elements included in the packet
extern int MakePermSprMatrixSendPacket (SparseMatrix spr, int *src_prm, int size, 
																					ptr_PacketNode pcknode);

// Prepare the structures required to receive some rows of a sparse matrix 
// and its corresponding permutation value
// * spr_aux refers to the SparseMatrix where the data will be stored
// * dst_prm includes the permutation vector
// * size is the number of used elements in dst_prm
// * pcknode, where the resulted packet appears
extern void MakePermSprMatrixRecvPacket (SparseMatrix spr_aux, int *dst_prm, int size, int weight, 
																					ptr_PacketNode pcknode);

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
extern int TestPacket (void *data);

// From the original data (spr,perm) and the elimination tree information (treetab,rangtab),
// the routine send the data related with a leaf to prc_dst, into the communicator comm.
// In lst appears the reference to the MPI_Isend operations which are not finalized. 
// The routine returns the number of nonzero elements which are sent
int SendLeafFromMatrix_old (SparseMatrix spr, int *perm, int *treetab, int *rangtab, 
													int prc_dst, int leaf, ptr_List lst, MPI_Comm comm);

// From the original matrix (spr) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst, into the communicator comm.
// In lst appears the reference to the MPI_Isend which are not finalized. 
// The routine returns the number of nonzero elements which are sent
int SendLeafFromMatrix (SparseMatrix spr, int prc_dst, int leaf, 
													ptr_ILU0Factor vFact, ptr_List lst, MPI_Comm comm);

// Compute the number of nonzeros elements of a PermSprMatrixRecvPacket packet
// * prc_src is the processor from which the messages is sent
// * sizes is the number of rows to be received
// * comm is the communicator in which the messages is sent
int ComputeSprMatrixWeights (int prc_src, int sizes, MPI_Comm comm);

// From elimination tree information (treetab,rangtab), the root sends
// the data related to a leaf:
// * sprL is the sparse matrix, and permL is the corresponding permutation.
// * headL is the vector which marks the begin and the end of each level (CSR way)
// * comm is the communicator in which the messages is received
// The routine returns the number of nonzero elements which are sent
extern int RecvLeafFromMatrix (int *treetab, int *rangtab, int root, int leaf,
													SparseMatrix sprL, int *permL, int *headL, MPI_Comm comm);

// The root sends the data from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
extern void RecvLeafFromMatrixF (ptr_ILU0Factor pFact, int root, int task, MPI_Comm comm);

// A slave sends the factor of task to root
// * comm is the communicator in which the messages is received
extern void SendFactorToMaster (ptr_ILU0Factor pFact, int task, int root, MPI_Comm comm);

// The root receives a factor from a slave (src)
// * vint includes the sizes of the packet 
// * comm is the communicator in which the messages is received
int RecvFactorFromSlaveFull (ptr_ILU0Factor vFact, int src, MPI_Comm comm);


// The root send a node to a slave (src)
// * comm is the communicator in which the messages is received
extern void SendNodeToSlave (ptr_ILU0Factor pFact, int src, int sizeTask, MPI_Comm comm);

// A slave receives a node from the root (root)
// * comm is the communicator in which the messages is received
extern void RecvNodeFromMaster (ptr_ILU0Factor pFact, int root, MPI_Comm comm);

/*********************************************************************************/

// Prepare the structures required to send some elements of 
// the original vector (vec) and its corresponding permutation value.
// * vec refers to the vector from where the data is obtained
// * src_prm includes the permutation vector
// * size is the number of used elements in src_prm
// * vvec is an auxiliar vector where the data is stored before to be sent
// * pcknode, where the resulted packet appears
extern void MakePermVectorSendPacket (double *vec, int *src_prm, int size, double *vvec, 
																			ptr_PacketNode pcknode);

// Prepare the structures required to receive some elements of a vector
// and its corresponding permutation value
// * vec_loc refers to the vector where the data is received
// * dst_prm includes the permutation vector
// * size is the number of used elements in src_prm
// * pcknode, where the resulted packet appears
extern void MakePermVectorRecvPacket (double *vec_loc, int *dst_prm, int size, 
																			ptr_PacketNode pcknode);

// From the original data (vector,perm) and the elimination tree information (treetab,rangtab),
// the routine send the data related with a leaf to prc_dst,
// where comm is the communicator in which the messages is sent
// In lst appears the reference to the MPI_Isend which are not finalized. 
extern void SendLeafFromVector (double *vec, int *perm, int *treetab, int *rangtab, 
													int prc_dst, int leaf, ptr_List lst, MPI_Comm comm);

// From the original vector (vec) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst,
// where comm is the communicator in which the messages is sent.
// In lst appears the reference to the MPI_Isend which are not finalized. 
extern void SendLeafFromVector2 (double *vec, int prc_dst, int leaf, 
													ptr_ILU0Factor vFact, ptr_List lst, MPI_Comm comm);

// From the original vector (vec) and the data included in *vFact,
// the routine send the data related with a leaf to prc_dst,
// where comm is the communicator in which the messages is sent.
// In lst appears the reference to the MPI_Isend which are not finalized. 
extern void SendLeafFromVector3 (double *vec, int prc_dst, int leaf, int vid,
													ptr_ILU0Factor pFact, ptr_List lst, MPI_Comm comm);

// From elimination tree information (treetab,rangtab), the root sends
// the data related to a leaf:
// * vecL is the vector, and permL is the corresponding permutation.
// * headL is the vector which marks the begin and the end of each level (CSR way)
// * comm is the communicator in which the messages is received
extern void RecvLeafFromVector (int *treetab, int *rangtab, int root, int leaf,
													double *vecL, int *permL, int *headL, MPI_Comm comm);

// From the original vector (vec), the root sends the corresponding elements 
// from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
extern void RecvLeafFromVectorF (ptr_ILU0Factor pFact, int root, int task, MPI_Comm comm);

// From the original vector (vec), the root sends the corresponding elements 
// from the task leaf, which are stored in pFact
// * comm is the communicator in which the messages is received
extern void RecvLeafFromVectorF3 (ptr_ILU0Factor pFact, int root, int vid, int task, 
																		int adjustVect, MPI_Comm comm);

void RecvLeafFromMatrixPCG (ptr_ILU0Factor pFact, int root, int task, MPI_Comm comm) ;
/*********************************************************************************/

extern void SendVectorResolutionTransform (ptr_ILU0Factor pFact, int prc_dst, int node, int UpDown,
																		int vidI1, int vidI2, int vidO2, ptr_List lst, MPI_Comm comm);

extern void RecvVectorResolutionTransform (ptr_ILU0Factor pFact, int root, int node, 
																		int vidI1, int vidI2, int vidO2, ptr_List lst, MPI_Comm comm);

/*********************************************************************************/

typedef struct {
	MPI_Comm comm_metis, comm_prec, comm_solver;
	int color, color2, root_m, root_p;
} Ilpck_Comm, *ptr_Ilpck_Comm;

//
extern void CreateILU0Communicator (ptr_Ilpck_Comm ptr_ilpk_comm, int root, int nprocs);

extern void RemoveILU0Communicator (ptr_Ilpck_Comm ptr_ilpk_comm);

/*********************************************************************************/

#endif
