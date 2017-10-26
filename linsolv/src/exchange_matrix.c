/*
 * exchange_matrix.c
 */

#include "exchange_matrix.h"

#include <mpi.h>
#include <stddef.h>
#include "util.h"
#include "setup_comm.h"

static int tag      = 47;
static int nrequest = 0;

static MPI_Request *request = NULL;
static MPI_Status  *status  = NULL;
static MPI_Comm comm        = MPI_COMM_WORLD;

static double **sendbuf = NULL;
static double **recvbuf = NULL;

/*******************************************************************************
*
*******************************************************************************/
void matrix_comm_init(CommMap *comap, int ncols)
{
  CHECK(comap != NULL);

  nrequest = comap->ncommdomains;

  request = check_malloc(2 * nrequest * sizeof(MPI_Request));
  status  = check_malloc(2 * nrequest * sizeof(MPI_Status));

  sendbuf = check_malloc(nrequest * sizeof(double*));
  recvbuf = check_malloc(nrequest * sizeof(double*));

  for(int i = 0; i < nrequest; i++)
  {
    const int k = comap->commpartner[i];
    const int sendcount = comap->sendcount[k];
    const int recvcount = comap->recvcount[k];

    if(sendcount > 0)
      sendbuf[i] = check_malloc(ncols * sendcount * sizeof(double));
    else
      sendbuf[i] = NULL;

    if(recvcount > 0)
      recvbuf[i] = check_malloc(ncols * recvcount * sizeof(double));
    else
      recvbuf[i] = NULL;
  }
} /* matrix_comm_init() */

/*******************************************************************************
*
*******************************************************************************/
void matrix_comm_end(void)
{
  check_free(request);
  request = NULL;
  check_free(status);
  status = NULL;

  for(int i = 0; i < nrequest; i++)
  {
    check_free(sendbuf[i]);
    check_free(recvbuf[i]);
  }

  check_free(sendbuf);
  sendbuf = NULL;
  check_free(recvbuf);
  recvbuf = NULL;

  nrequest = 0;
} /* matrix_comm_end() */

/*******************************************************************************
*
*******************************************************************************/
void exchange_matrix(const CommMap *comap, Matrix matrix)
{
  if(comap == NULL || comap->ndomains <= 1)
    return;

  const int ncols = matrix.cols;
  double **mat = matrix.m;
  const int ncommpartner = comap->ncommdomains;
  int *commpartner = comap->commpartner;
  int i, j, col;

  for(i = 0; i < ncommpartner; i++)
  {
    const int source = commpartner[i];
    const int recvcount = comap->recvcount[source];

    if(recvcount > 0)
      MPI_Irecv(recvbuf[i], ncols * recvcount, MPI_DOUBLE, source, tag, comm,
                request + i);
  }

  for(i = 0; i < ncommpartner; i++)
  {
    const int dest = commpartner[i];
    const int sendcount = comap->sendcount[dest];

    if(sendcount > 0)
    {
      int *sendindex = comap->sendindex[dest];

      /* copy data to sendbuffer */
      for(j = 0; j < sendcount; j++)
        for(col = 0; col < ncols; col++)
          sendbuf[i][j * ncols + col] = mat[sendindex[j]][col];

      MPI_Isend(sendbuf[i], ncols * sendcount, MPI_DOUBLE, dest, tag, comm,
                request + nrequest + i);
    }
  }

  MPI_Waitall(2 * nrequest, request, status);

  for(i = 0; i < ncommpartner; i++)
  {
    const int source = commpartner[i];
    const int recvcount = comap->recvcount[source];

    if(recvcount > 0)
    {
      int *recvindex = comap->recvindex[source];

      /* copy data from recvbuffer */
      for(j = 0; j < recvcount; j++)
        for(col = 0; col < ncols; col++)
          mat[recvindex[j]][col] = recvbuf[i][j * ncols + col];
    }
  }
} /* exchange_matrix() */
