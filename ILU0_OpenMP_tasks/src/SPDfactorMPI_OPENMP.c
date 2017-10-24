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
#include "TaskQueue.h"
#include "ToolsILU0.h"
#include "ILU0Factor.h"

#include "Lists.h"
#include "ToolsMPI.h"
#include "ToolsMPI_OPENMP.h"
#include "TreeOrderings.h"
#include "SPDfactorMPI_OPENMP.h"


/*********************************************************************************/

//#define VERBOSE 1
//#define CHANGE_ORDER 1
// #define SORT_BY_TIME 1
// #define LOAD_TASK_QUEUE 1
// #define WRITE_TASK_QUEUE 1

/*********************************************************************************/

// This routine computes an elimination tree with nleaves nodes of matrix spr. 
// After it creates the vector on which the preconditioner will be,
// filling it with the data related with the tree structure.
// At the end, the filled vector is returned
// The parameter index indicates if 0-indexing or 1-indexing is used.
// The parameter ilpkcomms includes the information related to the communicators 
// on which the computations will be made.
ptr_ILU0Factor InitILU0FactorizationMPIOPENMP (SparseMatrix spr, int index, int nleaves,
																								Ilpck_Comm ilpkcomms, char *matrix) {
	// Definition of the local variables
	int i, j, dim = spr.dim1, dimL, dimT = nleaves*2;
	int *perm = NULL, *rangtab = NULL, *treetab = NULL;
	int *sizetab = NULL, *chldtab = NULL, *nmchtab = NULL, *brthtab = NULL, *wgthtab = NULL;
	int *hgthtab = NULL, *sorttab = NULL, *ownrtab = NULL, *v_ord=NULL;
	matInts mTab;
	matDoubles tGlb;
	ILU0Factor *vFact = NULL;
	double tt, tt1, tt2, tge1, tge2, tgu1, tgu2;
	int ind =0, start, salto, stop, l;
	char permutation[30] = "", rang[30] = "", tree[30] = "", str_leaves[30];
	char *token;

	// Definition of the variables related to MPI 
	int my_id, numprocs, root; 
	MPI_Comm comm;
	
	// Definition of the variables numprocs and my_id
	comm = ilpkcomms.comm_metis; root = ilpkcomms.root_m;
	MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);

	// Creation of the shared structures
	if (my_id == root) CreateInts (&perm, dim); 
	dimT = 2*nleaves; //By definition
	CreateMatrixInts (&mTab, SIZE_MAT_TAB, dimT); InitInts (mTab[0], dimT*SIZE_MAT_TAB, 0, 0);
	CreateMatrixDoubles (&tGlb, 2, SIZE_TIM_GLB); InitDoubles (tGlb[0], 2*SIZE_TIM_GLB, 0.0, 0.0);

	rangtab = mTab[RANG]; InitInts (rangtab, dimT, -4, 0);
	treetab = mTab[TREE]; InitInts (treetab, dimT, -6, 0);
  ownrtab = mTab[OWNR]; InitInts (ownrtab, dimT, -7, 0);
	v_ord = mTab[LEV];		InitInts (v_ord, dimT, -1, 0);

	/* get the first token */
  token = strtok(matrix, "/");
   
  /* walk through other tokens */
  while( token != NULL ) 
  {
		if (token[0]=='A')
  		strcpy(matrix, token);
    token = strtok(NULL, "/");
  }

  sprintf(str_leaves, "_%d_", nleaves);
	strcat(matrix, str_leaves);
	strcat(permutation, matrix); strcat(permutation, "Perm.txt");
	strcat(rang, matrix);				 strcat(rang, "Rang.txt");
	strcat(tree, matrix);	       strcat(tree, "Tree.txt");

	MPI_Barrier (comm); reloj (&tge1, &tgu1); tt1 = MPI_Wtime ();
	if (my_id == root) {
		int i_aux, *aux = NULL;
		char file_aux[100];
		// Read the permutation
		sprintf(file_aux, "../Data/%s", permutation);
		i_aux = ReadInts (file_aux, &aux);
		if (i_aux != dim) {
			printf ("Error reading perm (%d,%d)\n", dim, i_aux);
			exit (-1);
		}
		CopyInts (aux, perm, dim); RemoveInts (&aux);
		// Read the size of the nodes
		sprintf(file_aux, "../Data/%s", rang);
		i_aux = ReadInts (file_aux, &aux);
		if (i_aux != dimT) {
			printf ("Error reading rangtab (%d,%d)\n", dimT, i_aux);
			exit (-1);
		}
		CopyInts (aux, rangtab, dimT); RemoveInts (&aux);
		// Read the elimination tree
		sprintf(file_aux, "../Data/%s", tree);
		i_aux = ReadInts (file_aux, &aux);
		if (i_aux != dimT) {
			printf ("Error reading treetab (%d,%d)\n", dimT, i_aux);
			exit (-1);
		}
		CopyInts (aux, treetab, dimT); RemoveInts (&aux);
	}
	reloj (&tge2, &tgu2); tt2 = MPI_Wtime (); tt = tt2 - tt1;
	tGlb[0][TGPMET] = tge2-tge1; tGlb[1][TGPMET] = tgu2-tgu1;

	MPI_Bcast (rangtab, dimT, MPI_INT, root, comm);
	MPI_Bcast (treetab, dimT, MPI_INT, root, comm);

	if (my_id == root) {
		printf ("Ordering Time (%2d) = %20.15e = (%20.15e,%20.15e)\n", my_id, tt, tge2-tge1, tgu2-tgu1);
#ifdef VERBOSE
		printf ("%s --> ", "rangtab"); PrintFInts (rangtab, dimT, 10, 0);
		printf ("%s --> ", "treetab"); PrintFInts (treetab, dimT, 10, 0);
#endif
	}

	if (ilpkcomms.color) {
		// Initialization of auxiliar vectors (mark and ownr are already initialized to 0)
		chldtab = mTab[CHLD]; nmchtab = mTab[NMCH]; brthtab = mTab[BRTH]; wgthtab = mTab[WGTH]; 
		ownrtab = mTab[OWNR]; hgthtab = mTab[HGTH]; sizetab = mTab[SIZE]; sorttab = mTab[SORT];
		// Compute sizetab, chldtab, nmchtab, brthtab, hgthtab vectors and the variable dimL (the number of nodes)
		dimL = ComputeEliminationTreeVectors (treetab, chldtab, nmchtab, brthtab, hgthtab, dimT);
		ComputeLengthfromHeader (rangtab, sizetab, dimL); 
#ifdef VERBOSE
		if (my_id == root) {
			printf ("header	 --> "); PrintFSeqInts	(0, dimL, 10, 0);
			printf ("sizetab --> "); PrintFInts (sizetab, dimL, 10, 0);
			printf ("hgthtab --> "); PrintFInts (hgthtab, dimL, 10, 0);
  		printf ("chldtab --> "); PrintFInts (chldtab, dimL, 10, 0);
      printf ("brthtab --> "); PrintFInts (brthtab, dimL, 10, 0);
		}
#endif
		// Compute the maximum number of nonzeros of each node (wgthtab vector)
		if (my_id == root) {
			ComputeWeightNodes (spr, index, perm, rangtab, sizetab, wgthtab, dimL);
		}
		MPI_Bcast (wgthtab, dimT, MPI_INT, root, comm);
#ifdef VERBOSE
		if (my_id == root) {
			printf ("wgthtab --> "); PrintFInts (wgthtab, dimL, 10, 0);
		}
#endif
		CopyInts (wgthtab, sorttab, dimL);
		for(i = 0; i < dimL; i++){
  		if(chldtab[i]==-1){
     		v_ord[ind]=i;
    		ind++;
   		}
		}
		start = 0;
		salto = nleaves;
		stop = nleaves;
		l = nleaves;
 		while ( l < dimL){
    	for(i = start; i < stop; i++){
    	  for( j = 0; j < dimL; j++){
    	    if(v_ord[i]==chldtab[j]){
    	      v_ord[ind] = j;
    	      ind++;
    	    }
    	  }
    	}
    	salto = salto/2;
    	start = stop;
    	stop = stop + salto;
    	l = l + salto;
  	}	
		// Create and initialize the vector of preconditioners
		vFact = CreateILU0FactorVector (dimL);
		vFact->indP = index; vFact->perm = perm ;
		vFact->mTab = mTab ; vFact->tGlb = tGlb ;
		vFact->dimL = dimL ; vFact->dimT = dimT ;
		for (i=1; i<dimL; i++) vFact[i] = *vFact;
	}
	// Return the vector of preconditioners
	return vFact;
}

#define NEW_MPIVERSION 1

// The routine computes the multilevel ILU factorization of the sparse matrix spr.
// The initial permutation is included in the files whose names appear as parameters.
int ILU0FactorizationMPIOPENMP (SparseMatrix spr, int index, int nleaves, paramFactor parFac, 
					Ilpck_Comm ilpkcomms, ptr_ILU0Factor *ppFact, char *matrix) { 
	 // Definition of the global vectors and variables
  int i;  
  int *treetab = NULL;
  int *chldtab = NULL, *brthtab = NULL;
  int remaining_tsks, *ownrtab = NULL, *hgthtab = NULL, *v_ord = NULL;
  ptr_ILU0Factor vFact;
  // Definition of the local vectors and variables
  matDoubles tGlb = NULL;
  matInts mTab;
  int shift;
  int my_id, numprocs, nleavespr;
  int indexM = 0, ierr = 0; 
  MPI_Comm comm;
  double tge1, tge2, tgu1, tgu2, tgt1, tgt2; // tt, tut1, tut2;

  // Initialization of the MPI variables
  srand (0);
  comm = ilpkcomms.comm_prec; 
  MPI_Comm_size(comm, &numprocs); MPI_Comm_rank(comm, &my_id);
  nleavespr = (nleaves / numprocs);
  if ((nleaves % numprocs) > 0) {
    printf ("Incorrect parameters in ILU0LeavesDistribution (%d,%d)\n", nleaves, numprocs);
    PrintTrace (); exit (-1);
  }

  // Computation of the partitioning
  MPI_Barrier (comm); reloj (&tge1, &tgu1); tgt1 = MPI_Wtime ();
  vFact = InitILU0FactorizationMPIOPENMP (spr, index, nleaves, ilpkcomms, matrix);
  mTab = vFact->mTab; tGlb = vFact->tGlb;
  reloj (&tge2, &tgu2); tgt2 = MPI_Wtime ();
  tGlb[0][TGFMET] = tge2-tge1; tGlb[1][TGFMET] = tgu2-tgu1;

  // Computation of the preconditioner
  tGlb[0][TGFPRC] =     -tge1; tGlb[1][TGFPRC] =     -tgu1;
  MPI_Barrier(comm); reloj (&tge1, &tgu1); tgt1 = MPI_Wtime ();
  if (ilpkcomms.color) {
    // Definition of the local variables of each thread
    int tid = omp_get_thread_num(), tsk, chld;
    int dst, src, flag = 1, cond, proc;
    double te = 0.0, tle1, tle2, tlu1, tlu2;
    int OMP_MPI = 1;
    int **dependencies, child, n_deps=3, father;
    // Initialization of auxiliar vectors (mark and ownr are already initialized to 0)
    treetab = mTab[TREE];
    chldtab = mTab[CHLD]; brthtab = mTab[BRTH]; 
    ownrtab = mTab[OWNR]; hgthtab = mTab[HGTH]; 
    v_ord = vFact->mTab[LEV];
#ifdef VERBOSE
    if (my_id == 0) printf ("(BEGIN) ILU0LeavesDistribution\n");
#endif
    ILU0LeavesDistribution_OPENMP (spr, index, nleaves, parFac, ilpkcomms, vFact);
#ifdef VERBOSE
    if (my_id == 0) printf (" (END)  ILU0LeavesDistribution\n");
#endif
    remaining_tsks = 2 * nleavespr - 1;
    dependencies = (int**) malloc (remaining_tsks * sizeof(int*));
    for (tsk=0; tsk<remaining_tsks; tsk++){
  	dependencies[tsk] = (int*) malloc (n_deps * sizeof(int));
  	child=chldtab[tsk];
  	for(i=0; i < n_deps-1; i++){
    		dependencies[tsk][i]=(child==-1)?tsk:child;
    		if(child!=-1)
    			child=brthtab[child];
    	}
    	father = treetab[tsk];
    	dependencies[tsk][n_deps-1]=(father==-1)?tsk:father;
    }
    // Initialize the number of processes to be processed
   	MPI_Barrier(comm); reloj (&tge1, &tgu1); tgt1 = MPI_Wtime ();
   	shift = v_ord[my_id*nleavespr];
		#pragma omp parallel
  	{	
   	#pragma omp single
   	{
   	for (i = 0; i < remaining_tsks; i++){
		#pragma omp task depend(in:chldtab[dependencies[i][0]], chldtab[dependencies[i][1]]) depend(out:chldtab[i]) firstprivate(i) private(tsk)
		{
    	double te = 0.0, tle1, tle2, tlu1, tlu2;
      // Try to get an actived and non-processed node of the tree
			tsk = i + shift;
			tid = omp_get_thread_num();
      reloj (&tle1, &tlu1); //tut1 = MPI_Wtime ();
      if (vFact[tsk].tLoc == NULL) {
      	CreateMatrixDoubles (&vFact[tsk].tLoc, 2, SIZE_TIM_LOC);
        InitDoubles (vFact[tsk].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
      }		
      // Compute the matrix to be factorized
      ownrtab[tsk] = tid; chld = chldtab[i];
      if (chld != -1) { // The new node is not a leaf
      	te = BuildNodeFromChildren (indexM, vFact, tsk);
      }
      reloj (&tle2, &tlu2); //tut2 = MPI_Wtime ();
      vFact[tsk].tLoc[0][TLACUM] += te;
      vFact[tsk].tLoc[0][TLBLDP] += (tle2-tle1) - te;
      vFact[tsk].tLoc[1][TLBLDP] += tlu2-tlu1;
      reloj (&tle1, &tlu1);
     	// Computation of the multilevel ILU preconditioner
      FactorILU0SparseMatrix (vFact+tsk, parFac, tsk, (chld == -1));
      if (vFact[tsk].ierr == 0) {
      	vFact[tsk].dimX = vFact[tsk].sprM.dim1;
       	vFact[tsk].dimF = vFact[tsk].sprF.dim1;
       	reloj (&tle2, &tlu2);
       	vFact[tsk].tLoc[0][TLFCTP] += tle2-tle1;
       	vFact[tsk].tLoc[1][TLFCTP] += tlu2-tlu1;
			}
		}
  	}
  	#pragma omp taskwait
  	}
 		}
  	MPI_Barrier(comm);
  	tsk = (i-1) + shift; 
  	remaining_tsks = 0;			
  	// MPI code
  	if ((remaining_tsks == 0) && (OMP_MPI)) {
  		OMP_MPI = 0; remaining_tsks = hgthtab[tsk] - 1;
  		flag = 1; cond = 1; proc = my_id;
  	}
  	if ((remaining_tsks > 0) && (OMP_MPI == 0)){
			dst = proc |  flag; src = proc & ~flag;
      cond = (proc & flag) > 0; flag <<= 1;
      if (proc == src) {
				SendFactorToMaster (vFact+tsk, tsk, dst, comm);
        cond = 0; remaining_tsks = 0;
			} else {
        reloj (&tle1, &tlu1); //tut1 = MPI_Wtime ();
        tsk = RecvFactorFromSlaveFull (vFact, src, comm);
				if (vFact[tsk].tLoc == NULL) {
		    	CreateMatrixDoubles (&vFact[tsk].tLoc, 2, SIZE_TIM_LOC);
        	InitDoubles (vFact[tsk].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
        }
        reloj (&tle2, &tlu2); //tut2 = MPI_Wtime ();
        vFact[tsk].tLoc[0][TLACUM] += 0.0;
        vFact[tsk].tLoc[0][TLBLDP] += (tle2-tle1);
        vFact[tsk].tLoc[1][TLBLDP] += (tlu2-tlu1);
			}	
		}
  	while ((remaining_tsks > 0)) { 
  		tsk = treetab[tsk];
			reloj (&tle1, &tlu1); //tut1 = MPI_Wtime ();
			if (vFact[tsk].tLoc == NULL) {
	    	CreateMatrixDoubles (&vFact[tsk].tLoc, 2, SIZE_TIM_LOC);
        InitDoubles (vFact[tsk].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
      }
      // Compute the matrix to be factorized
      ownrtab[tsk] = tid; chld = chldtab[tsk];
			if (chld != -1) {
				te = BuildNodeFromChildren(indexM, vFact, tsk); // Build the node from the factors computed for the children of node tsk
			}
			reloj (&tle2, &tlu2); //tut2 = MPI_Wtime ();
      vFact[tsk].tLoc[0][TLACUM] += te;
      vFact[tsk].tLoc[0][TLBLDP] += (tle2-tle1) - te;
      vFact[tsk].tLoc[1][TLBLDP] += tlu2-tlu1;

      reloj (&tle1, &tlu1); //tut1 = MPI_Wtime ();
			//Factorize
			FactorILU0SparseMatrix (vFact+tsk, parFac, tsk, (chld == -1)); 
      reloj (&tle2, &tlu2); //tut2 = MPI_Wtime ();
			vFact[tsk].dimX = vFact[tsk].sprM.dim1; vFact[tsk].dimF = vFact[tsk].sprF.dim1;
			vFact[tsk].tLoc[0][TLFCTP] += tle2-tle1;
      vFact[tsk].tLoc[1][TLFCTP] += tlu2-tlu1;

			remaining_tsks--;
			
	 		if ((remaining_tsks > 0) && (OMP_MPI == 0)) {
				dst = proc |  flag; src = proc & ~flag;
        cond = (proc & flag) > 0; flag <<= 1;
        if (proc == src) {
					SendFactorToMaster (vFact+tsk, tsk, dst, comm);
        	cond = 0; remaining_tsks = 0;
				} else{
		    	reloj (&tle1, &tlu1); //tut1 = MPI_Wtime ();
      		tsk = RecvFactorFromSlaveFull (vFact, src, comm);
			 		if (vFact[tsk].tLoc == NULL) {
		      	CreateMatrixDoubles (&vFact[tsk].tLoc, 2, SIZE_TIM_LOC);
          	InitDoubles (vFact[tsk].tLoc[0], 2*SIZE_TIM_LOC, 0.0, 0.0);
        	}
        	reloj (&tle2, &tlu2); //tut2 = MPI_Wtime ();
        	vFact[tsk].tLoc[0][TLACUM] += 0.0;
        	vFact[tsk].tLoc[0][TLBLDP] += (tle2-tle1);
        	vFact[tsk].tLoc[1][TLBLDP] += (tlu2-tlu1);
				}
			}
		}	
  }
 
  // Print the time of the preconditioner computation
  reloj (&tge2, &tgu2); tgt2 = MPI_Wtime ();
  tGlb[0][TGFPRC] += tge2     ; tGlb[1][TGFPRC] += tgu2     ;
  tGlb[0][TGPPRC]  = tge2-tge1; tGlb[1][TGPPRC]  = tgu2-tgu1;
	if(my_id == 0){
  	printf ("Preconditioner Time (%d)= (%20.15e,%20.15e,%20.15e)\n", my_id, tge2-tge1, tgu2-tgu1, tgt2-tgt1);
	}
#ifdef VERBOSE
  // Print Preconditioner Time
  printf ("Preconditioner Time = (%20.15e,%20.15e)\n", tge2-tge1, tgu2-tgu1);
  // Print ownrtab
  printf ("ownrtab --> "); PrintFInts (vFact->mTab[OWNR], vFact->dimL, 10, 0);
#endif
  // Return data
  *ppFact = vFact;

  return ierr;
}

