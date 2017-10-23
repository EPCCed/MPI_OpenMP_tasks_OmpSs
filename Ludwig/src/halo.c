
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pe.h"
#include "coords.h"
#include "model.h"
#include "lb_model_s.h"
#include "targetDP.h"
#include "io_harness.h"
#include "control.h"

/*****************************************************************************
 *
 *  halo_planes
 *
 *  Sends 6 MPI messages (the planes) for the Non-Blocking version.
 *
 *****************************************************************************/

void halo_planes_tasks(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo, indexreal;
  int count;
  int nlocal[3];

  const int tagf = 900;
  const int tagb = 901;

  /* The ranks of neighbouring planes */
  int pforwX, pbackX, pforwY, pbackY, pforwZ, pbackZ;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  /* Allocate size of sendplanes and number of elements send with each plane */
  int nsendXZ, nsendYZ, nsendXY;
  nsendYZ = NVEL*lb->ndist*nlocal[Y]*nlocal[Z];
  nsendXZ = NVEL*lb->ndist*nlocal[X]*nlocal[Z];
  nsendXY = NVEL*lb->ndist*nlocal[X]*nlocal[Y];

  /* Allocate message sizes for plane send/receives */
  lb->hl .sendforwYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .sendbackYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .recvforwYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .recvbackYZ = (double *) malloc(nsendYZ*sizeof(double));
  lb->hl .sendforwXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .sendbackXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .recvforwXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .recvbackXZ = (double *) malloc(nsendXZ*sizeof(double));
  lb->hl .sendforwXY = (double *) malloc(nsendXY*sizeof(double));
  lb->hl .sendbackXY = (double *) malloc(nsendXY*sizeof(double));
  lb->hl .recvforwXY = (double *) malloc(nsendXY*sizeof(double));
  lb->hl .recvbackXY = (double *) malloc(nsendXY*sizeof(double));

  /* Receive planes in the X-direction */
  /* PPM, NMM, P=Positive, M=Middle, N=Negative for the XYZ directions respectively */
  pforwX = nonblocking_cart_neighb(PMM);
  pbackX = nonblocking_cart_neighb(NMM);

  MPI_Irecv(&lb->hl.recvforwYZ[0], nsendYZ, MPI_DOUBLE, pforwX, tagb, comm,
            &lb->hl.recvreq[0]);
  MPI_Irecv(&lb->hl .recvbackYZ[0], nsendYZ, MPI_DOUBLE, pbackX, tagf, comm,
            &lb->hl.recvreq[1]);

  /* Receive planes in the Y-direction */
  pforwY = nonblocking_cart_neighb(MPM);
  pbackY = nonblocking_cart_neighb(MNM);

  MPI_Irecv(&lb->hl .recvforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagb, comm,
            &lb->hl.recvreq[2]);
  MPI_Irecv(&lb->hl .recvbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagf, comm,
            &lb->hl.recvreq[3]);

  /* Receive planes in the Z-direction */
  pforwZ = nonblocking_cart_neighb(MMP);
  pbackZ = nonblocking_cart_neighb(MMN);

  MPI_Irecv(&lb->hl.recvforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagb, comm,
            &lb->hl.recvreq[4]);
  MPI_Irecv(&lb->hl.recvbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagf, comm,
            &lb->hl.recvreq[5]);

  /* Send in the X-direction (YZ plane) */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(nlocal[X], jc, kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl.sendforwYZ[count] = fptr[indexreal];

          index = coords_index(1, jc, kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl.sendbackYZ[count] = fptr[indexreal];
          ++count;
        }
      }
    }
  }
  assert(count == nsendYZ);

  MPI_Issend(&lb->hl .sendbackYZ [0], nsendYZ, MPI_DOUBLE, pbackX, tagb, comm,
             &lb->hl.sendreq[0]);
  MPI_Issend(&lb->hl .sendforwYZ [0], nsendYZ, MPI_DOUBLE, pforwX, tagf, comm,
             &lb->hl.sendreq[1]);


  /* Send in the Y-direction (XZ plane) */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (kc = 1; kc <= nlocal[Z]; kc++) {

          index = coords_index(ic, nlocal[Y], kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendforwXZ[count] = fptr[indexreal];

          index = coords_index(ic, 1, kc);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendbackXZ[count] = fptr[indexreal];
          ++count;
        }
      }
    }
  }
  assert(count == nsendXZ);

  MPI_Issend(&lb->hl .sendbackXZ[0], nsendXZ, MPI_DOUBLE, pbackY, tagb, comm,
             &lb->hl.sendreq[2]);
  MPI_Issend(&lb->hl .sendforwXZ[0], nsendXZ, MPI_DOUBLE, pforwY, tagf, comm,
             &lb->hl.sendreq[3]);


  /* Finally, Send in the Z-direction (XY plane) */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {
        for (jc = 1; jc <= nlocal[Y]; jc++) {

          index = coords_index(ic, jc, nlocal[Z]);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendforwXY[count] = fptr[indexreal];

          index = coords_index(ic, jc, 1);
          indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
          lb->hl .sendbackXY[count] = fptr[indexreal];
          ++count;
        }
      }
    }
  }
  assert(count == nsendXY);

  MPI_Issend(&lb->hl .sendbackXY[0], nsendXY, MPI_DOUBLE, pbackZ, tagb, comm,
             &lb->hl.sendreq[4]);
  MPI_Issend(&lb->hl .sendforwXY[0], nsendXY, MPI_DOUBLE, pforwZ, tagf, comm,
             &lb->hl.sendreq[5]);

  return;

}

/*****************************************************************************
 *
 *  halo_edges
 *
 *  Sends 12 MPI messages (the edges) for the Non-Blocking version.
 *
 *****************************************************************************/


void halo_edges_tasks(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int index, indexhalo, indexreal;
  int count;
  int nlocal[3];

  const int tagnn = 903;
  const int tagnp = 904;
  const int tagpn = 905;
  const int tagpp = 906;

  /* Ranks of neighbouring edges.
     Xpp, X is the direction, parallel to which edges are
     sent pp refers to the YZ directions respectively*/
  int Xpp, Xpn, Xnp, Xnn, Ypp, Ypn, Ynp, Ynn, Zpp, Zpn, Znp, Znn;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr; /*pointer to the distribution*/
  if (get_step()) /* we are in the timestep, so use target copy */
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;  /* we are not in a timestep, use host copy */

  int nsendX, nsendY, nsendZ;
  nsendX = NVEL*lb->ndist*nlocal[X];
  nsendY = NVEL*lb->ndist*nlocal[Y];
  nsendZ = NVEL*lb->ndist*nlocal[Z];

  /* Allocate message sizes for edges send/receives */
  lb->hl .sendXnn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXnn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendXnp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXnp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendXpn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXpn = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendXpp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .recvXpp = (double *) malloc(nsendX*sizeof(double));
  lb->hl .sendYnn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYnn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .sendYnp = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYnp = (double *) malloc(nsendY*sizeof(double));

  lb->hl .sendYpn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYpn = (double *) malloc(nsendY*sizeof(double));
  lb->hl .sendYpp = (double *) malloc(nsendY*sizeof(double));
  lb->hl .recvYpp = (double *) malloc(nsendY*sizeof(double));
  lb->hl .sendZnn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZnn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .sendZnp = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZnp = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .sendZpn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZpn = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .sendZpp = (double *) malloc(nsendZ*sizeof(double));
  lb->hl .recvZpp = (double *) malloc(nsendZ*sizeof(double));

  /* Receive edges parallel to x-direction*/
  Xnn = nonblocking_cart_neighb(MNN);
  Xnp = nonblocking_cart_neighb(MNP);
  Xpn = nonblocking_cart_neighb(MPN);
  Xpp = nonblocking_cart_neighb(MPP);

  MPI_Irecv(&lb->hl.recvXnn[0], nsendX, MPI_DOUBLE, Xnn, tagpp, comm, &lb->hl.recvreq[6]);
  MPI_Irecv(&lb->hl.recvXnp[0], nsendX, MPI_DOUBLE, Xnp, tagpn, comm, &lb->hl.recvreq[7]);
  MPI_Irecv(&lb->hl.recvXpn[0], nsendX, MPI_DOUBLE, Xpn, tagnp, comm, &lb->hl.recvreq[8]);
  MPI_Irecv(&lb->hl.recvXpp[0], nsendX, MPI_DOUBLE, Xpp, tagnn, comm, &lb->hl.recvreq[9]);

  /* Receive edges parallel to y-direction*/
  Ynn = nonblocking_cart_neighb(NMN);
  Ynp = nonblocking_cart_neighb(NMP);
  Ypn = nonblocking_cart_neighb(PMN);
  Ypp = nonblocking_cart_neighb(PMP);

  MPI_Irecv(&lb->hl.recvYnn[0], nsendY, MPI_DOUBLE, Ynn, tagpp, comm, &lb->hl.recvreq[10]);
  MPI_Irecv(&lb->hl.recvYnp[0], nsendY, MPI_DOUBLE, Ynp, tagpn, comm, &lb->hl.recvreq[11]);
  MPI_Irecv(&lb->hl.recvYpn[0], nsendY, MPI_DOUBLE, Ypn, tagnp, comm, &lb->hl.recvreq[12]);
  MPI_Irecv(&lb->hl.recvYpp[0], nsendY, MPI_DOUBLE, Ypp, tagnn, comm, &lb->hl.recvreq[13]);

  /* Receive edges parallel to z-direction*/
  Znn = nonblocking_cart_neighb(NNM);
  Znp = nonblocking_cart_neighb(NPM);
  Zpn = nonblocking_cart_neighb(PNM);
  Zpp = nonblocking_cart_neighb(PPM);

  MPI_Irecv(&lb->hl.recvZnn[0], nsendZ, MPI_DOUBLE, Znn, tagpp, comm, &lb->hl.recvreq[14]);
  MPI_Irecv(&lb->hl.recvZnp[0], nsendZ, MPI_DOUBLE, Znp, tagpn, comm, &lb->hl.recvreq[15]);
  MPI_Irecv(&lb->hl.recvZpn[0], nsendZ, MPI_DOUBLE, Zpn, tagnp, comm, &lb->hl.recvreq[16]);
  MPI_Irecv(&lb->hl.recvZpp[0], nsendZ, MPI_DOUBLE, Zpp, tagnn, comm, &lb->hl.recvreq[17]);

  /* Send edges parallel to x-direction */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (ic = 1; ic <= nlocal[X]; ic++) {

        index = coords_index(ic, 1, 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXnn[count] = fptr[indexreal];

        index = coords_index(ic, 1, nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXnp[count] = fptr[indexreal];

        index = coords_index(ic, nlocal[Y], 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXpn[count] = fptr[indexreal];

        index = coords_index(ic, nlocal[Y], nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendXpp[count] = fptr[indexreal];
        ++count;
      }
    }
  }
  assert(count == nsendX);

  MPI_Issend(&lb->hl .sendXpp[0], nsendX, MPI_DOUBLE, Xpp, tagpp, comm, &lb->hl.sendreq[6]);
  MPI_Issend(&lb->hl .sendXpn[0], nsendX, MPI_DOUBLE, Xpn, tagpn, comm, &lb->hl.sendreq[7]);
  MPI_Issend(&lb->hl .sendXnp[0], nsendX, MPI_DOUBLE, Xnp, tagnp, comm, &lb->hl.sendreq[8]);
  MPI_Issend(&lb->hl .sendXnn[0], nsendX, MPI_DOUBLE, Xnn, tagnn, comm, &lb->hl.sendreq[9]);

  /* Send edges parallel to y-direction */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (jc = 1; jc <= nlocal[Y]; jc++) {

        index = coords_index(1, jc, 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYnn[count] = fptr[indexreal];

        index = coords_index(1, jc, nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYnp[count] = fptr[indexreal];

        index = coords_index(nlocal[X], jc, 1);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYpn[count] = fptr[indexreal];

        index = coords_index(nlocal[X], jc, nlocal[Z]);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendYpp[count] = fptr[indexreal];
        ++count;
      }
    }
  }
  assert(count == nsendY);

  MPI_Issend(&lb->hl .sendYpp[0], nsendY, MPI_DOUBLE, Ypp, tagpp, comm, &lb->hl.sendreq[10]);
  MPI_Issend(&lb->hl .sendYpn[0], nsendY, MPI_DOUBLE, Ypn, tagpn, comm, &lb->hl.sendreq[11]);
  MPI_Issend(&lb->hl .sendYnp[0], nsendY, MPI_DOUBLE, Ynp, tagnp, comm, &lb->hl.sendreq[12]);
  MPI_Issend(&lb->hl .sendYnn[0], nsendY, MPI_DOUBLE, Ynn, tagnn, comm, &lb->hl.sendreq[13]);

  /* Send edges parallel to z-direction */
  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {
      for (kc = 1; kc <= nlocal[Z]; kc++) {

        index = coords_index(1, 1, kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZnn[count] = fptr[indexreal];

        index = coords_index(1, nlocal[Y], kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZnp[count] = fptr[indexreal];

        index = coords_index(nlocal[X], 1, kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZpn[count] = fptr[indexreal];

        index = coords_index(nlocal[X], nlocal[Y], kc);
        indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
        lb->hl .sendZpp[count] = fptr[indexreal];
        ++count;
      }
    }
  }
  assert(count == nsendZ);

  MPI_Issend(&lb->hl .sendZpp[0] , nsendZ, MPI_DOUBLE, Zpp, tagpp, comm, &lb->hl.sendreq[14]);
  MPI_Issend(&lb->hl .sendZpn[0], nsendZ, MPI_DOUBLE, Zpn, tagpn, comm, &lb->hl.sendreq[15]);
  MPI_Issend(&lb->hl .sendZnp[0], nsendZ, MPI_DOUBLE, Znp, tagnp, comm, &lb->hl.sendreq[16]);
  MPI_Issend(&lb->hl .sendZnn[0], nsendZ, MPI_DOUBLE, Znn, tagnn, comm, &lb->hl.sendreq[17]);

  return;
}

/*****************************************************************************
 *
 *  halo_corners
 *
 *  Sends 8 MPI messages (the corners) for the Non-Blocking version.
 *
 *****************************************************************************/

void halo_corners_tasks(lb_t * lb) {

  int ic, jc, kc;
  int n, p;
  int count;
  int index, indexhalo, indexreal;
  int nlocal[3];

  const int tagnnn = 907;
  const int tagnnp = 908;
  const int tagnpn = 909;
  const int tagnpp = 910;
  const int tagpnn = 911;
  const int tagpnp = 912;
  const int tagppn = 913;
  const int tagppp = 914;

  /* Ranks of neighbouring corners XYZ direction*/
  int ppp, ppn, pnp, pnn, npp, npn, nnp, nnn;

  MPI_Comm comm = cart_comm();

  assert(lb);

  coords_nlocal(nlocal);

  double* fptr;
  if (get_step())
    fptr=lb->tcopy->f;
  else
    fptr=lb->f;

  int nsend;
  nsend = NVEL*lb->ndist*1;

  /* Allocate message sizes for plane send/receives */
  lb->hl .sendnnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendnnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendnpn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendnpp = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendpnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendpnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendppn = (double *) malloc(nsend*sizeof(double));
  lb->hl .sendppp = (double *) malloc(nsend*sizeof(double));

  lb->hl .recvnnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvnnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvnpn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvnpp = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvpnn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvpnp = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvppn = (double *) malloc(nsend*sizeof(double));
  lb->hl .recvppp = (double *) malloc(nsend*sizeof(double));

  nnn = nonblocking_cart_neighb(NNN);
  nnp = nonblocking_cart_neighb(NNP);
  npn = nonblocking_cart_neighb(NPN);
  npp = nonblocking_cart_neighb(NPP);
  pnn = nonblocking_cart_neighb(PNN);
  pnp = nonblocking_cart_neighb(PNP);
  ppn = nonblocking_cart_neighb(PPN);
  ppp = nonblocking_cart_neighb(PPP);

  MPI_Irecv(&lb->hl.recvnnn[0], nsend, MPI_DOUBLE, nnn, tagppp, comm, &lb->hl.recvreq[18]);
  MPI_Irecv(&lb->hl.recvnnp[0], nsend, MPI_DOUBLE, nnp, tagppn, comm, &lb->hl.recvreq[19]);
  MPI_Irecv(&lb->hl.recvnpn[0], nsend, MPI_DOUBLE, npn, tagpnp, comm, &lb->hl.recvreq[20]);
  MPI_Irecv(&lb->hl.recvnpp[0], nsend, MPI_DOUBLE, npp, tagpnn, comm, &lb->hl.recvreq[21]);
  MPI_Irecv(&lb->hl.recvpnn[0], nsend, MPI_DOUBLE, pnn, tagnpp, comm, &lb->hl.recvreq[22]);
  MPI_Irecv(&lb->hl.recvpnp[0], nsend, MPI_DOUBLE, pnp, tagnpn, comm, &lb->hl.recvreq[23]);
  MPI_Irecv(&lb->hl.recvppn[0], nsend, MPI_DOUBLE, ppn, tagnnp, comm, &lb->hl.recvreq[24]);
  MPI_Irecv(&lb->hl.recvppp[0], nsend, MPI_DOUBLE, ppp, tagnnn, comm, &lb->hl.recvreq[25]);

  count = 0;
  for (p = 0; p < NVEL; p++) {
    for (n = 0; n < lb->ndist; n++) {

      index = coords_index(1, 1, 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnnn[count] = fptr[indexreal];

      index = coords_index(1, 1, nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnnp[count] = fptr[indexreal];

      index = coords_index(1, nlocal[Y], 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnpn[count] = fptr[indexreal];

      index = coords_index(1, nlocal[Y], nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendnpp[count] = fptr[indexreal];

      index = coords_index(nlocal[X], 1, 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendpnn[count] = fptr[indexreal];

      index = coords_index(nlocal[X], 1, nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendpnp[count] = fptr[indexreal];

      index = coords_index(nlocal[X], nlocal[Y], 1);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendppn[count] = fptr[indexreal];

      index = coords_index(nlocal[X], nlocal[Y], nlocal[Z]);
      indexreal = LB_ADDR(lb->nsite, lb->ndist, NVEL, index, n, p);
      lb->hl .sendppp[count] = fptr[indexreal];
      count++;
    }
  }

  MPI_Issend(&lb->hl .sendppp[0], nsend, MPI_DOUBLE, ppp, tagppp, comm, &lb->hl.sendreq[18]);
  MPI_Issend(&lb->hl .sendppn[0], nsend, MPI_DOUBLE, ppn, tagppn, comm, &lb->hl.sendreq[19]);
  MPI_Issend(&lb->hl .sendpnp[0], nsend, MPI_DOUBLE, pnp, tagpnp, comm, &lb->hl.sendreq[20]);
  MPI_Issend(&lb->hl .sendpnn[0], nsend, MPI_DOUBLE, pnn, tagpnn, comm, &lb->hl.sendreq[21]);
  MPI_Issend(&lb->hl .sendnpp[0], nsend, MPI_DOUBLE, npp, tagnpp, comm, &lb->hl.sendreq[22]);
  MPI_Issend(&lb->hl .sendnpn[0], nsend, MPI_DOUBLE, npn, tagnpn, comm, &lb->hl.sendreq[23]);
  MPI_Issend(&lb->hl .sendnnp[0], nsend, MPI_DOUBLE, nnp, tagnnp, comm, &lb->hl.sendreq[24]);
  MPI_Issend(&lb->hl .sendnnn[0], nsend, MPI_DOUBLE, nnn, tagnnn, comm, &lb->hl.sendreq[25]);

  return;
}
