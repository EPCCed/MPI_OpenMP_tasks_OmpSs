/*
 * test_linsolv.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include "util.h"
#include "setup_comm.h"
#include "matrix.h"
#include "read_data.h"
#include "linsys.h"
#include "linsolv.h"
#include "exchange_matrix.h"

#ifdef SET_NSWEEPS
# define N_SWEEPS SET_NSWEEPS
#else
# define N_SWEEPS 5
#endif

int main(int argc, char *argv[])
{
  int num_sweeps = N_SWEEPS;
  int nsolver = 4;
  BSLinSysMethod solverIDs[4] = {BSLS_POINT_IMPLICIT,
                                 BSLS_JACOBI,
                                 BSLS_GAUSS_SEIDEL,
                                 BSLS_SYMM_GAUSS_SEIDEL};

  char filename[1024];
  int myrank, nprocs;

  /* init mpi and comap */
  mpi_parallel_init(argc, argv, &myrank, &nprocs);
  CommMap *comap = init_comap(myrank, nprocs);

  int nthreads = omp_get_max_threads();

  if(DBG_PRINT || myrank == 0)
    printf("[%d] Number of threads: %d\n"
           "[%d] Number of procs:   %d\n",
           myrank, nthreads, myrank, nprocs);

  if(!getenv("OMP_NUM_THREADS"))
    printf("[%d] OMP_NUM_THREADS is not set!\n", myrank);

  if(DBG_PRINT)
  {
    printf("[%d] INFO: dynamic threads are %s\n",
           myrank, omp_get_dynamic() ? "enabled" : "disabled");
    printf("[%d] INFO: nested parallel regions are %s\n",
           myrank, omp_get_nested() ? "enabled" : "disabled");

#   pragma omp parallel default(none) shared(myrank, nprocs, nthreads)
    {
      int threadid = omp_get_thread_num();

      printf("[%d] Hello from thread: %3d of %3d in rank %4d of %4d (MPI)\n",
             myrank, threadid, nthreads, myrank, nprocs);
    }
  }

  /* read and calc communication data only in parallel case */
  if(nprocs > 1)
  {
    /* read comm data */
    if(DBG_PRINT || myrank == 0)
      printf("[%d] Read communication data\n", myrank);
    sprintf(filename, "comap_domain_%dof%d", myrank, nprocs);
    read_comap(filename, comap);

    /* compute comm tables */
    create_recvsend_index(comap);
  } /* parallel case */

  /* read matrix */
  BlockSparseMatrix bsmatrix;
  if(DBG_PRINT || myrank == 0)
    printf("[%d] Read block sparse matrix\n", myrank);
  sprintf(filename, "linsys_bsmatrix_domain_%dof%d", myrank, nprocs);
  read_block_sparse_matrix(filename, &bsmatrix);

  /* read rhs */
  Matrix rhs;
  if(DBG_PRINT || myrank == 0)
    printf("[%d] Read rhs\n", myrank);
  sprintf(filename, "linsys_rhs_domain_%dof%d", myrank, nprocs);
  read_matrix(filename, &rhs);

  /* check dimensions */
  if(nprocs > 1)
  {
    CHECK(comap->nownpoints == bsmatrix.num_rows);
    CHECK(comap->nownpoints + comap->naddpoints == rhs.rows);
  }
  else
    CHECK(bsmatrix.num_rows == rhs.rows);
  CHECK(bsmatrix.block_size_row == rhs.cols);
  CHECK(bsmatrix.block_size_row == bsmatrix.block_size_row);

  /* allocate mpi buffers, requests and statuses */
  if(nprocs > 1)
    matrix_comm_init(comap, rhs.cols);

  /* setup linear system */
  BSMLinSys bsm_linsys = {0};
  bsm_linsys.A = bsmatrix;
  bsm_linsys.rhs = rhs;
  bsm_linsys.swap = check_malloc(bsmatrix.num_rows * sizeof(int*));
  for(int row = 0; row < bsmatrix.num_rows; row++)
    bsm_linsys.swap[row] = check_calloc(bsmatrix.block_size_row - 1,
                                        sizeof(int));
  bsm_linsys.decomposition_state = NOT_DECOMPOSED;
  bsm_linsys.diag_first = check_diag_entry_is_first(&bsmatrix);

  if(!bsm_linsys.diag_first)
    printf("[%d] Block sparse matrix NOT in expected format!\n", myrank);

  Matrix soln = generateMatrix(rhs.rows, rhs.cols);

  /* print task version infos */
  if(DBG_PRINT || myrank == 0)
  {
    printf("[%d] Active task version: %s\n", myrank, TASKVERSION);
  }
  if(DBG_PRINT && USE_TASKLOOP == 0)
  {
    printf("[%d]  used chunksize: %d\n", myrank, SET_CHUNKSIZE);
  }

  if(DBG_PRINT || myrank == 0)
    printf("[%d] Number of sweeps: %d\n", myrank, num_sweeps);

  /* decompose matrix */
  if(DBG_PRINT || myrank == 0)
    printf("[%d] Decompose matrix\n", myrank);
  linsys_lu_decomp_diag(&bsm_linsys);

  /* linear solver loop */
  for(int i = 0; i < nsolver; i++)
  {
    BSLinSysMethod lin_solver = solverIDs[i];

    if(DBG_PRINT || myrank == 0)
      printf("[%d] Test solver '%s'\n", myrank, get_linsolv_name(lin_solver));

    linsolv(lin_solver, num_sweeps, &bsm_linsys, soln, comap);

    print_colum_checksums_6(soln, bsmatrix.num_rows);

    clearMatrix(soln); /* reset */
  } /* linear solver loop */

  /* free test data */
  delete_bsm_linear_system(bsm_linsys); /* frees also bsmatrix and rhs */
  deleteMatrix(soln);

  if(DBG_PRINT || myrank == 0)
    printf("[%d] All done\n", myrank);

  /* free comm ressources */
  free_comap(&comap);
  matrix_comm_end();
  mpi_parallel_end();
} /* main() */
