/*
 * setup_comm.c
 */

#include "setup_comm.h"

#include <stddef.h> /* NULL */
#include <mpi.h>
#include "util.h"

/*******************************************************************************
*
*******************************************************************************/
CommMap *init_comap(int thisdom, int ndom)
{
  CommMap *ptr = (CommMap *)check_malloc(sizeof(CommMap));

  ptr->thisdomain = thisdom;
  ptr->ndomains = ndom;
  ptr->ncommdomains = 0;
  ptr->nownpoints = 0;
  ptr->naddpoints = 0;

  ptr->addpoint_owner = NULL;
  ptr->addpoint_id = NULL;
  ptr->commpartner = NULL;
  ptr->sendcount = NULL;
  ptr->recvcount = NULL;
  ptr->recvindex = NULL;
  ptr->sendindex = NULL;

  return ptr;
} /* init_comap() */

/*******************************************************************************
*
*******************************************************************************/
void free_comap(CommMap **comap)
{
  CommMap *old = *comap;

  if(old == NULL)
    return;

  check_free(old->addpoint_owner);
  check_free(old->addpoint_id);
  check_free(old->commpartner);
  check_free(old->sendcount);
  check_free(old->recvcount);

  if(old->sendindex != NULL)
    for(int i = 0; i < old->ndomains; i++)
      check_free(old->sendindex[i]);
  check_free(old->sendindex);

  if(old->recvindex != NULL)
    for(int i = 0; i < old->ndomains; i++)
      check_free(old->recvindex[i]);
  check_free(old->recvindex);

  check_free(*comap);
} /* free_comap() */

/*******************************************************************************
*
*******************************************************************************/
void create_recvsend_index(CommMap *comap)
{
  int tag = 3001;
  int nown, nadd, ndom, myid, i;

  if(comap == NULL || comap->ncommdomains == 0 || comap->ndomains <= 1)
    return;

  CHECK(comap->commpartner != NULL);
  CHECK(comap->sendcount != NULL);
  CHECK(comap->recvcount != NULL);

  if(comap->naddpoints > 0)
  {
    CHECK(comap->addpoint_owner != NULL);
    CHECK(comap->addpoint_id != NULL);
  }

  /*----------------------------------------------------------------------------
  | init
  ----------------------------------------------------------------------------*/
  nown = comap->nownpoints;
  nadd = comap->naddpoints;
  ndom = comap->ndomains;
  myid = comap->thisdomain;

  if(comap->sendindex != NULL)
    for(i = 0; i < ndom; i++)
      check_free(comap->sendindex[i]);
  check_free(comap->sendindex);

  if(comap->recvindex != NULL)
    for(i = 0; i < ndom; i++)
      check_free(comap->recvindex[i]);
  check_free(comap->recvindex);

  comap->sendindex = check_malloc(ndom * sizeof(int *));
  comap->recvindex = check_malloc(ndom * sizeof(int *));

  for(i = 0; i < ndom; i++)
    comap->sendindex[i] = comap->recvindex[i] = NULL;

  for(i = 0; i < comap->ncommdomains; i++)
  {
    int j;
    int k = comap->commpartner[i];

    comap->sendindex[k] = check_malloc(comap->sendcount[k] * sizeof(int));
    for(j = 0; j < comap->sendcount[k]; j++)
      comap->sendindex[k][j] = -1;

    if(comap->recvcount[k] > 0)
    {
      int count = 0;

      comap->recvindex[k] = check_malloc(comap->recvcount[k] * sizeof(int));
      for(j = 0; j < nadd; j++)
        if(comap->addpoint_owner[j] == k)
          comap->recvindex[k][count++] = nown + j;
    }
  }

  /*----------------------------------------------------------------------------
  | buffer
  ----------------------------------------------------------------------------*/
  size_t sz = 0;
  for(i = 0; i < comap->ncommdomains; i++)
  {
    int k = comap->commpartner[i];
    sz = MAX(sz, (size_t)comap->sendcount[k]);
    sz = MAX(sz, (size_t)comap->recvcount[k]);
  }

  int *ibuf = check_malloc(sz * sizeof(int));

  /*----------------------------------------------------------------------------
  | exchange
  ----------------------------------------------------------------------------*/
  for(i = 0; i < comap->ncommdomains; i++)
  {
    int k          = comap->commpartner[i];
    int recvcount  = comap->recvcount[k];
    int sendcount  = comap->sendcount[k];
    int *recvindex = comap->recvindex[k];
    int *sendindex = comap->sendindex[k];
    int j;

    if(k > myid) /* first send */
    {
      for(j = 0; j < recvcount; j++)
        ibuf[j] = comap->addpoint_id[recvindex[j] - nown];

      MPI_Send(ibuf, recvcount * sizeof(int), MPI_BYTE, k, tag, MPI_COMM_WORLD);

      MPI_Recv(ibuf, sendcount * sizeof(int), MPI_BYTE, k, tag, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      for(j = 0; j < sendcount; j++)
        sendindex[j] = ibuf[j];
    }
    else  /* first receive */
    {
      MPI_Recv(ibuf, sendcount * sizeof(int), MPI_BYTE, k, tag, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      for(j = 0; j < sendcount; j++)
        sendindex[j] = ibuf[j];
      for(j = 0; j < recvcount; j++)
        ibuf[j] = comap->addpoint_id[recvindex[j] - nown];

      MPI_Send(ibuf, recvcount * sizeof(int), MPI_BYTE, k, tag, MPI_COMM_WORLD);
    }
  } /* for(i = 0; i < comap->ncommdomains; i++) */

  check_free(ibuf);
} /* create_recvsend_index() */

/*******************************************************************************
*
*******************************************************************************/
void mpi_parallel_init(int argc, char *argv[], int *rank, int *ndom)
{
  int provided, required;

#ifdef USE_MPI_MULTI_THREADED
  required = MPI_THREAD_MULTIPLE;
#else
  required = MPI_THREAD_SERIALIZED;
#endif

  MPI_Init_thread(&argc, &argv, required, &provided);
  CHECK(provided >= required);

  MPI_Comm_rank(MPI_COMM_WORLD, rank);
  MPI_Comm_size(MPI_COMM_WORLD, ndom);

  DBG_MSG("[%d] Hello from rank %4d of %4d (MPI)\n", *rank, *rank, *ndom);
} /* mpi_parallel_init() */

/*******************************************************************************
*
*******************************************************************************/
void mpi_parallel_end(void)
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
} /* mpi_parallel_end() */
