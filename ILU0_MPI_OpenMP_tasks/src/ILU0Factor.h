#ifndef ILU0FactorType

#define ILU0FactorType 1

#include <ScalarVectors.h>
#include <SparseMatricesNew.h>

/**********************************************************************/

// Definition of the datatype used to store the parameters which are required for the computations
typedef struct paramFactor {
	int nparts, nthreads, maxit;
} paramFactor, *ptr_paramFactor;

// Number of the parameters of the executions
#define NUMBER_PARAMS 7

// This routine verify the number of parameters and its contents
extern void VerificationParameters (int argc, char **argv, int *nleaves, paramFactor *parFac);

/**********************************************************************/

// Definition of the datatype for the preconditioner
typedef struct factorType {
	SparseMatrix sprLDU, UF;           // The LDU factorization and the rectangular matrices r    elated to ILU0
	int ierr;                          // The result of the computation of the preconditioner

	int indP, indM, indF, indR;				 // Define the indexing (0 or 1) of perm, sprM, sprF and sprR.
	int *perm;                         // Parallel ND
	matInts mTab;											 // Set of vectors to Parallel ND
	int dimL, dimT;                    // Sizes of the rangtab and treetab vectors
	int dimX, dimF;                    // Sizes of original sparse matrix and of the computed Schur complement
	int nlev, *headL;                  // Size of the path, and CSR way vector from begin to end
	double *divs;                      // Divisor from a leaf to the root
	SparseMatrix sprM, sprF;           // Original matrix and factor
	int dimM;                          // Size of the matrix 
	int dimV, tskV;                    // Addition of the sizes of the leaves and task where the malloc is made
	matDoubles mDia;                   // Set of vectors on which the diagonal vectors of the leading blocks are
	int *permM;                        // For the leaves, permutation vector to get/put data
	SparseMatrix sprR;								 // Matrix on which the PCG process is applied
	matDoubles mPCG;									 // Set of vectors which are used to compute the PCG
	matDoubles mVcL;									 // Set of vectors which are used to compute the resolution and the transformation
	int *ntsks;                        // Vector of number of tasks in each thread
	matDoubles tGlb, tLoc;						 // Set of vectors where the computational times are stored
} ILU0Factor, *ptr_ILU0Factor;

/**********************************************************************/

// Constants related to the global counters
#define SIZE_TIM_GLB  6
#define TGFMET        0
#define TGPMET        1
#define TGFPRC        2
#define TGPPRC        3
#define TGFPCG        4
#define TGPPCG        5

// Constants related to the local counters
#define SIZE_TIM_LOC 17
#define TLBLDP        0
#define TLACUM        1
#define TLFCTP        2
#define TLAMGP        3

#define TLBLDU        4
#define TLRESU        5
#define TLSOLU        6

#define TLBLDD        7
#define TLRESD        8
#define TLSOLD        9

#define TLPROD       10
#define TLVEC0       11
#define TLVEC1       12
#define TLVEC2       13
#define TLVEC3       14
#define TLVEC4       15
#define TLVEC5       16

/**********************************************************************/

// Constants related to the vectors which describe 
// the structure of the elimination tree
#define SIZE_MAT_TAB 14
#define RANG          0
#define TREE          1
#define SIZE          2
#define CHLD          3
#define NMCH          4
#define BRTH          5
#define WGTH          6
#define MARK          7
#define OWNR          8
#define HGTH          9
#define SORT         10
#define SORU				 11	
#define SORD				 12
#define LEV					 13	

// Constants used to manage the diagonals
#define SIZE_MAT_DIA  3
#define DIAG          0
#define DIA1          1
#define DIA2          2

// Constants used to apply the preconditioner
#define SIZE_MAT_VCL  8
#define VECL          0
#define BUFL          1
#define SOLL          2
#define AUXL          3
#define AUXL0         3
#define AUXL1         4
#define AUXL2         5
#define TRNL          6

// Constants used to the compute the PCG
#define SIZE_MAT_PCG  8
#define VECT_BU       0
#define VECT_XC       1
#define VECT_ZU       2
#define VECT_RU       3
#define VECT_YC       4
#define VECT_DC       5
#define VECT_RC       6
#define VECT_MATD     7
#define VECT_VECL    10
#define VECT_BUFL    11
#define VECT_SOLL    12
#define VECT_AUXL    13
#define VECT_AUXL0   13
#define VECT_AUXL1   14
#define VECT_AUXL2   15
#define VECT_TRNF    16
#define VECT_SCAL    20

/*********************************************************************************/

// Create a vector of preconditioners with dim elements
extern ptr_ILU0Factor CreateILU0FactorVector (int dim);

// Clean a vector of preconditioners, maintaining only the preconditioner,
// the permutation and the METIS structures.
void CleanILU0FactorVector (ptr_ILU0Factor *ppFact);

// Remove a vector of preconditioners 
extern void RemoveILU0FactorVector (ptr_ILU0Factor *ppFact);

// Write a vector of preconditioners 
extern void WriteILU0FactorVector (ptr_ILU0Factor vFact, char *nameMatr, 
																						char *nameFact, char *nameScal);
// This routine returns the accumulation of the sizes of the leaves: 
// - The addition of the number of rows
// - The addition of the number of nonzeros
extern void GetSizesLeaves (ptr_ILU0Factor vFact, int *dimML, int *nnzML);

// This routine returns the accumulation of the sizes of the leaves: 
// - The addition of the number of rows
// - The addition of the number of nonzeros
extern void GetSizesLeavesPCG (ptr_ILU0Factor vFact, int *dimRL, int *nnzRL);

/*********************************************************************************/

// This routine returns the address of the vector vid related to the structure pFact
extern double *IdentifyVectorResolution (ptr_ILU0Factor pFact, int vid);

/****************************************************************************/

// This routine creates the vector tRes, whose size is equal the the number of leaves, 
// and copies the contents of a local timer to this new vector.
extern void GetLocalTimer (ptr_ILU0Factor vFact, int timer, double **tRes);

// This routine prints the information included in the counters delimited by
// nTim1 and nTim2, acording to the value of typeP:
// * < 1, it prints all the information, from the nodes to the final summarize.
// * < 2, it prints the information, from the threads to the final summarize.
// * < 3, it only prints the summarize.
// First, the information is accumulated by each level of the threads on the 
// matrix mSumThr, and after the addition and the maximum of each level is computed.
// Finally the addition of the maximum is computed to obtain the estimation
// of the parallel cost, and the final addition to obtain the sequential cost.
// These additions accumulate the counters by subsets (n1,n2) and (n3,n4),
// and the results are returned in the parameters auxS and auxP.
extern void PrintBlockTime (int typeP, ptr_ILU0Factor vFact, int nTim1, int nTim2,
														matDoubles mSumMax, int hmax, matDoubles *mSumThr, int nthr,
														int n1, int n2, int n3, int n4, double *auxS, double *auxP);

// This routine prints the information of the counters, acording to the value of typeP:
// * < 1, it prints all the information, from the nodes to the final summarize.
// * < 2, it prints the information, from the threads to the final summarize.
// * < 3, it only prints the summarize.
extern void PrintTimesILU0Factor (int typeP, ptr_ILU0Factor vFact, int itrEnd);

/****************************************************************************/

#endif

